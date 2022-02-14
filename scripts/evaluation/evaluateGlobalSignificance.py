#!/usr/bin/env python
import argparse
import re
import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.optimize import root_scalar, minimize
from statsmodels.distributions.empirical_distribution import ECDF


parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("--method", choices=["asymptotics", "toys"], default="asymptotics")
parser.add_argument("--q0-sampling-distributions", nargs="*")
args = parser.parse_args()


import ROOT as R
R.gROOT.SetBatch(True)
R.gROOT.SetStyle("ATLAS")


# Load toy results
if args.method == "toys":
    assert args.q0_sampling_distributions is not None \
        and len(args.q0_sampling_distributions) == 20

    q0_toys = {}
    q0_ecdf = {}
    for fn in args.q0_sampling_distributions:
        m = re.search(r"q0_(\d+)\.csv$", fn)
        if not m:
            raise RuntimeError(f"Cannot parse mass from: {fn}")

        mass, = m.groups()
        mass = int(mass)

        q0_toys[mass] = pd.read_csv(fn)
        q0_ecdf[mass] = ECDF(q0_toys[mass]["q0"])


# Configuration
obs_z0 = 3.012651697087006
obs_pval = 1.0 - stats.norm.cdf(obs_z0)

if args.method == "toys":
    obs_q0 = obs_z0**2
    obs_pval = 1 - q0_ecdf[1000](obs_q0)
    obs_z0 = stats.norm.ppf(q0_ecdf[1000](obs_q0))

# Read global significance toys
df = pd.read_csv(args.infile, low_memory=False)

# Convert to proper dtypes
df = df.astype({
    "uncond_status": "int64",
    "cond_status": "int64",
    "uncond_covQual": "int64",
    "cond_covQual": "int64",
})

# Rename 'index' to 'toyindex' to avoid confusion with the index of
# the dataframe
if "index" in df.columns:
    df.rename(columns={"index": "toyindex"}, inplace=True)

# Check that there are no duplicates
assert not df.duplicated(["toyindex", "mass"]).any()

# Number of toys with unique index
total_toys = len(np.unique(df["toyindex"]))
print(f"Total number of toys: {total_toys}")

# Failed fits
df["failed_fit"] = (df["uncond_status"] != 0) | (df["cond_status"] != 0)
num_failed = df["failed_fit"].sum()
frac_failed = num_failed / len(df)
print(f"Failed fits: {num_failed} ({100 * frac_failed:.1f} %)")

# Transform q0 to one-sided discovery test statistic
df.loc[df["muhat"] <= 0, "q0"] = 0.0

# Only keep toys where all fits succeeded
df["good_toy"] = ~df.groupby("toyindex")["failed_fit"].transform(np.any)
print("Fraction of good toys: {:.1f} %".format(
    100 * df["good_toy"].sum() / len(df)))

df_good = df.loc[df["good_toy"]].copy()

# Set slightly negative q0 to 0 and check that none left
threshold = -0.05
df_good.loc[(df_good["q0"] < 0) & (df_good["q0"] > threshold), "q0"] = 0.0
#TODO: Remove me
df_good.loc[df_good["q0"] < 0, "q0"] = 0.0

assert len(df_good.loc[df_good["q0"] < 0]) == 0

# Calculate *local* significance
df_good["sig_asymptotics"] = np.sqrt(df_good["q0"])
df_good["sig_toys"] = 0.0

if args.method == "toys":
    for mass in df_good["mass"].unique():
        q0 = df_good.loc[df_good["mass"] == mass, "q0"]
        sig = stats.norm.ppf(q0_ecdf[mass](q0))

        # Getting infinities if q0 > q0 value seen for local toys
        mask_inf = np.isinf(sig)
        num_inf = np.count_nonzero(mask_inf)

        if num_inf > 0:
            print(f"Encountered infinities for mass: {mass}")
            print("Setting significance to 4.99")
            # This is fine as long as the default value > the observed
            # threshold which is ~3
            sig[mask_inf] = 4.99

        assert not np.any(np.isinf(sig))
        df_good.loc[df_good["mass"] == mass, "sig_toys"] = sig


# Maximum local significance
df_zmax = df_good.groupby("toyindex")["sig_asymptotics"].max().to_frame("max_sig_asymptotics")
df_zmax["max_sig_toys"] = df_good.groupby("toyindex")["sig_toys"].max()

# To calculate 1sigma interval of p-value
def binom_mle_interval(k, n, bracket=(0.01, 0.1)):
    llh = lambda p, k, n: k * np.log(p) + (n - k) * np.log(1 - p)
    delta_llh = lambda p: llh(p, k, n) - llh(k / n, k, n)

    root_lo = root_scalar(lambda p: -2 * delta_llh(p) - 1, bracket=[bracket[0], k / n])
    root_hi = root_scalar(lambda p: -2 * delta_llh(p) - 1, bracket=[k / n, bracket[1]])
    assert root_lo.converged and root_hi.converged

    return root_lo.root, root_hi.root

# Calculate global p-value / significance
num_exceeding = (df_zmax[f"max_sig_{args.method}"] > obs_z0).sum()
num_toys = len(df_zmax)
global_pval = num_exceeding / num_toys
global_sig = stats.norm.ppf(1 - global_pval)

print(f"Global p-value ({args.method}): {100 * global_pval:.2f} %")
print(f"Global significance ({args.method}): {global_sig:.2f}")

# Calculate confidence intervals on global p and sig
global_pval_lo, global_pval_hi = binom_mle_interval(num_exceeding, num_toys)
print(f"68% CL interval on p-global: [{100 * global_pval_lo:.2f} %, {100 * global_pval_hi:.2f} %]")

global_sig_lo, global_sig_hi = stats.norm.ppf(1 - global_pval_hi), stats.norm.ppf(1 - global_pval_lo)
print(f"68% CL interval on sig-global: [{global_sig_lo:.2f}, {global_sig_hi:.2f}]")


def fit_trial_factor(sig_max):
    # Constants removed
    log_pdf = lambda x, n: np.log(n) + (n - 1) * stats.norm.logcdf(x) - x**2 / 2
    nll = lambda n: -log_pdf(sig_max, n).sum()

    res = minimize(nll, x0=[15.], bounds=[[10., 21.]])
    if not res.success:
        print(res)

    assert res.success

    trial_factor, = res.x

    return trial_factor


# Fit the trial factor
trial_factor = fit_trial_factor(df_zmax[f"max_sig_{args.method}"])
tf_global_sig = stats.norm.ppf(1 - (1 - (1 - obs_pval)**trial_factor))
print(f"Trial factor: {trial_factor:.1f}")
print(f"Global significance (trial factor): {tf_global_sig:.2f}")

# PDF for trial factor for plot
tf_pdf = lambda x, n: n * stats.norm.cdf(x)**(n - 1) * np.exp(-x**2/2.) / np.sqrt(2 * np.pi)
x_tf = np.linspace(0, 5, 200)
norm = (5 / 100.) * len(df_zmax) # Bin width * number of toys
y_tf = norm * tf_pdf(x_tf, trial_factor)

g_tf = R.TGraph(len(x_tf), x_tf, y_tf)




# Plot
h_zmax = R.TH1F("h_zmax", "", 100, 0, 5)
for index, zmax in df_zmax[f"max_sig_{args.method}"].iteritems():
    h_zmax.Fill(zmax)

h_zmax.GetXaxis().SetTitle("Maximum local significance [#sigma]")
h_zmax.GetYaxis().SetTitle("Toy experiments")
h_zmax.SetMinimum(0)
h_zmax.SetMaximum(1.3 * h_zmax.GetMaximum())
h_zmax.SetLineColor(R.kBlue)
h_zmax.SetMarkerSize(0)

g_tf.SetLineWidth(2)
g_tf.SetLineColor(R.kOrange + 1)

latex = R.TLatex()
latex.SetNDC()
latex.SetTextFont(43)
latex.SetTextSize(19)

c = R.TCanvas("c", "", 800, 600)
h_zmax.Draw("HIST,E0")
g_tf.Draw("C,SAME")

# Labels
latex.DrawLatex(0.54, 0.85, f"Number of toys: {len(df_zmax)}")
latex.DrawLatex(0.54, 0.80, "Toys with Z_{0}^{max} > Z_{0}^{max,obs}: " + f"{num_exceeding}")
latex.DrawLatex(0.54, 0.72,
                "Global p-value (toys): {:.4f}#splitline{{#plus {:.4f}}}{{#minus {:.4f}}}".format(
                    global_pval, # Central value
                    global_pval_hi - global_pval, # up error
                    global_pval - global_pval_lo)) # down error
latex.DrawLatex(0.54, 0.64,
                "Global significance (toys): {:.2f}#splitline{{#plus {:.2f}}}{{#minus {:.2f}}}".format(
                    global_sig, # Central value
                    global_sig_hi - global_sig, # up error
                    global_sig - global_sig_lo)) # down error
latex.DrawLatex(0.54, 0.56, f"Trial factor: {trial_factor:.1f} (Z_{{global}} = {tf_global_sig:.2f})")

# Line at observed significance
line = R.TLine(obs_z0, 0, obs_z0, h_zmax.GetMaximum())
line.SetLineStyle(R.kDashed)
line.Draw()
latex.SetNDC(0)

if args.method == "toys":
    latex.DrawLatex(obs_z0 + 0.05, 0.2 * h_zmax.GetMaximum(),
                    "Z_{{0}}^{{max,obs}}(toys) = {:.2f}".format(obs_z0))
else:
    latex.DrawLatex(obs_z0 + 0.05, 0.2 * h_zmax.GetMaximum(),
                    "Z_{{0}}^{{max,obs}}(asymptotics) = {:.2f}".format(obs_z0))

c.RedrawAxis()
c.SaveAs("zmax.pdf")
