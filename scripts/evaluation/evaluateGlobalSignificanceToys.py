#!/usr/bin/env python
import argparse
import re
import warnings

from joblib import Parallel, delayed
from scipy import stats
import numpy as np
import pandas as pd

warnings.simplefilter(action='ignore', category=FutureWarning)
import statsmodels.api as sm

import ROOT as R
R.gROOT.SetBatch(True)
R.gROOT.SetStyle("ATLAS")

from common import load_toys, binom_mle_interval, fit_trial_factor
from common import ToyPvalueCalculator


parser = argparse.ArgumentParser()
parser.add_argument("toys")
parser.add_argument("--q0-sampling-distributions", nargs="+", required=True)
args = parser.parse_args()


# Load toy experiments (global significance)
df = load_toys(args.toys)

total_toys = len(np.unique(df["toyindex"]))
total_failed = df["failed_fit"].sum()

df_good = df.loc[df["good_toy"]].copy()
total_good = len(np.unique(df_good["toyindex"]))
df_good.loc[df_good["q0"] < 0, "q0"] = 0.0

print(f"Total number of toys: {total_toys}")
print(f"Failed fits: {total_failed}")
print(f"Good toys: {total_good}")

num_neg_q0 = np.count_nonzero(df_good["q0"] < 0)
if num_neg_q0 > 0:
    print(f"Warning: Encountered {num_neg_q0} fits with negative q0")


# Load q0 sampling distributions (null)
toy_pval_calc = ToyPvalueCalculator()
assert len(args.q0_sampling_distributions) == 20

for fn in args.q0_sampling_distributions:
    m = re.search(r"q0_(\d+)\.csv$", fn)
    if not m:
        raise RuntimeError(f"Cannot parse mass from: {fn}")

    mass, = m.groups()
    mass = int(mass)

    q0 = pd.read_csv(fn)
    toy_pval_calc.add_q0_distribution(mass, q0)


# Observed results
obs_q0 = 3.012651697087006**2  # Significance with asymptotics squared
obs_pval = toy_pval_calc.get_pval(1000, obs_q0)
obs_z0 = stats.norm.ppf(1 - obs_pval)

print(f"Obs. q0: {obs_q0:.4f}")
print(f"Obs. p-value: {obs_pval:.4f}")
print(f"Obs. significance: {obs_z0:.4f}")


# Nominal result
df_good["pval"] = toy_pval_calc.get_pval(df_good["mass"], df_good["q0"])
df_good["sig"] = stats.norm.ppf(1 - df_good["pval"])

# Check for infinities
mask_inf = np.isinf(df_good["sig"])
num_inf = np.count_nonzero(mask_inf)
if num_inf > 0:
    print(f"Encountered infinities for mass: {mass} "
          "-- Setting to sqrt(q0)")
    df_good.loc[mask_inf, "sig"] = np.sqrt(df_good.loc[mask_inf, "q0"])


# Calculate global significance
df_zmax = df_good.groupby("toyindex")["sig"].max().to_frame("max_sig")

num_exceeding = (df_zmax["max_sig"] > obs_z0).sum()
num_toys = len(df_zmax)
global_pval = num_exceeding / num_toys
global_sig = stats.norm.ppf(1 - global_pval)

print(f"Global p-value: {global_pval:.4f}")
print(f"Global significance: {global_sig:.4f}")

global_pval_lo, global_pval_hi = binom_mle_interval(num_exceeding, num_toys)
print("68% CL interval on p-global: "
      f"[{100 * global_pval_lo:.2f} %, {100 * global_pval_hi:.2f} %]")

global_sig_lo, global_sig_hi = \
    stats.norm.ppf(1 - global_pval_hi), stats.norm.ppf(1 - global_pval_lo)
print(f"68% CL interval on sig-global: "
      f"[{global_sig_lo:.2f}, {global_sig_hi:.2f}]")


# Trial factor
trial_factor = fit_trial_factor(df_zmax["max_sig"])
tf_global_sig = stats.norm.ppf(1 - (1 - (1 - obs_pval)**trial_factor))
print(f"Trial factor: {trial_factor:.1f}")
print(f"Global significance (trial factor): {tf_global_sig:.2f}")


def tf_pdf(x, n):
    return (n * stats.norm.cdf(x)**(n - 1) * np.exp(-x**2/2.)
            / np.sqrt(2 * np.pi))


x_tf = np.linspace(0, 5, 200)
norm = (5 / 100.) * len(df_zmax)  # Bin width * number of toys
y_tf = norm * tf_pdf(x_tf, trial_factor)
g_tf = R.TGraph(len(x_tf), x_tf, y_tf)


# Plot
h_zmax = R.TH1F("h_zmax", "", 100, 0, 5)
h_zmax.FillN(len(df_zmax), df_zmax["max_sig"].values, R.nullptr)

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
# g_tf.Draw("C,SAME")

latex.DrawLatex(0.54, 0.85, f"Number of toys: {len(df_zmax)}")
latex.DrawLatex(0.54, 0.80, "Toys with Z_{0}^{max} > Z_{0}^{max,obs}: "
                f"{num_exceeding}")
latex.DrawLatex(0.54, 0.72,
                "Global p-value (toys): "
                "{:.4f}#splitline{{#plus {:.4f}}}{{#minus {:.4f}}}".format(
                    global_pval,  # Central value
                    global_pval_hi - global_pval,  # up error
                    global_pval - global_pval_lo))  # down error
latex.DrawLatex(0.54, 0.64,
                "Global significance (toys): "
                "{:.2f}#splitline{{#plus {:.2f}}}{{#minus {:.2f}}}".format(
                    global_sig,  # Central value
                    global_sig_hi - global_sig,  # up error
                    global_sig - global_sig_lo))  # down error
# latex.DrawLatex(0.54, 0.56, "Trial factor: "
#                 f"{trial_factor:.1f}")

# Line at observed significance
line = R.TLine(obs_z0, 0, obs_z0, h_zmax.GetMaximum())
line.SetLineStyle(R.kDashed)
line.Draw()
latex.SetNDC(0)

latex.DrawLatex(obs_z0 + 0.05, 0.2 * h_zmax.GetMaximum(),
                "Z_{{0}}^{{max,obs}}(toys) = {:.2f}".format(obs_z0))

c.RedrawAxis()
c.SaveAs("zmax_toys.pdf")
del c


# Bootstrapping
def bootstrap(rng, num_bootstraps):
    zglobal_bootstraps = []

    for i in range(num_bootstraps):
        # Resample q0 values in p value calculator
        toy_pval_calc_resampled = toy_pval_calc.bootstrap(rng)

        # Recalculate global significance

        # Observed results (needs to be recalculated since q0 toys changed)
        obs_pval = toy_pval_calc_resampled.get_pval(1000, obs_q0)
        obs_z0 = stats.norm.ppf(1 - obs_pval)

        df_good["pval"] = \
            toy_pval_calc_resampled.get_pval(df_good["mass"], df_good["q0"])
        df_good["sig"] = stats.norm.ppf(1 - df_good["pval"])

        # Check for infinities
        mask_inf = np.isinf(df_good["sig"])
        num_inf = np.count_nonzero(mask_inf)
        if num_inf > 0:
            df_good.loc[mask_inf, "sig"] = np.sqrt(df_good.loc[mask_inf, "q0"])

        # Calculate global significance
        df_zmax = df_good.groupby("toyindex")["sig"].max().to_frame("max_sig")

        num_exceeding = (df_zmax["max_sig"] > obs_z0).sum()
        num_toys = len(df_zmax)
        global_pval = num_exceeding / num_toys
        global_sig = stats.norm.ppf(1 - global_pval)

        zglobal_bootstraps.append(global_sig)

    return zglobal_bootstraps


seed_seq = np.random.SeedSequence(202596828575806468932732570400374383977)
child_seq = seed_seq.spawn(8)

zglobal_bootstraps = Parallel(n_jobs=8, verbose=10)(
    delayed(bootstrap)(np.random.default_rng(seq), 12500) for seq in child_seq)

zglobal_bootstraps = np.array(zglobal_bootstraps)
zglobal_bootstraps = zglobal_bootstraps.ravel()


bootstrap_mean = zglobal_bootstraps.mean()
bootstrap_std = zglobal_bootstraps.std(ddof=1)

print("Global significance (bootstrap mean): "
      f"{bootstrap_mean:.3f}")
print("Global significance error (bootstrap std): "
      f"{bootstrap_std:.3f}")
print("Error on mean global significance: "
      f"{bootstrap_std / np.sqrt(len(zglobal_bootstraps)):.3f}")


# Density estimate
zmax_kde = sm.nonparametric.KDEUnivariate(zglobal_bootstraps)
zmax_kde.fit()

x = np.linspace(1.85, 2.15, 100)
y = zmax_kde.evaluate(x)

g_kde = R.TGraph(len(x), x, y)
g_kde.SetLineWidth(2)
g_kde.SetLineColor(R.kBlue)

h_dummy = R.TH1F("h_dummy", "", 100, 1.85, 2.15)
h_dummy.GetXaxis().SetTitle(r"Z_{global}")
h_dummy.GetYaxis().SetTitle(r"Probability Density")
h_dummy.SetMinimum(0)
h_dummy.SetMaximum(1.3 * y.max())

latex = R.TLatex()
latex.SetNDC()
latex.SetTextFont(43)
latex.SetTextSize(19)

leg = R.TLegend(0.7, 0.8, 0.9, 0.9)
leg.SetTextFont(43)
leg.SetTextSize(19)
leg.SetBorderSize(0)
leg.AddEntry(g_kde, "KDE", "l")

c = R.TCanvas("c", "", 800, 600)
h_dummy.Draw("AXIS")
g_kde.Draw("L,SAME")

latex.DrawLatex(0.25, 0.8, f"Mean: {bootstrap_mean:.3f}")
latex.DrawLatex(0.25, 0.75, f"Std. Dev.: {bootstrap_std:.3f}")
leg.Draw()

c.RedrawAxis()
c.SaveAs("zmax_bootstrap.pdf")
c.SaveAs("zmax_bootstrap.root")
