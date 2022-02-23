#!/usr/bin/env python
import argparse

import numpy as np
import scipy.stats as stats

from common import load_toys, binom_mle_interval, fit_trial_factor


parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("--replace-failures", action="store_true",
                    help="Replace failing fits with q0 = 0")
args = parser.parse_args()


import ROOT as R
R.gROOT.SetBatch(True)
R.gROOT.SetStyle("ATLAS")


# Configuration
obs_z0 = 3.012651697087006
obs_pval = 1.0 - stats.norm.cdf(obs_z0)

df = load_toys(args.infile)


# Print out statistics
total_toys = len(np.unique(df["toyindex"]))
num_failed = df["failed_fit"].sum()
frac_failed = num_failed / len(df)
print(f"Total number of toys: {total_toys}")
print(f"Failed fits: {num_failed} ({100 * frac_failed:.1f} %)")
print("Fraction of good toys: {:.1f} %".format(
    100 * df["good_toy"].sum() / len(df)))

if args.replace_failures:
    print("Replacing failed fits with q0 = 0...")
    df.loc[df["failed_fit"], "q0"] = 0.0
    df["good_toy"] = True

df_good = df.loc[df["good_toy"]].copy()


# Set negative q0 to 0
df_good.loc[df_good["q0"] < 0, "q0"] = 0.0
assert len(df_good.loc[df_good["q0"] < 0]) == 0


# Calculate *local* significance
df_good["sig"] = np.sqrt(df_good["q0"])


# Maximum local significance
df_zmax = df_good.groupby("toyindex")["sig"].max().to_frame("max_sig")


# Calculate global p-value / significance
num_exceeding = (df_zmax["max_sig"] > obs_z0).sum()
num_toys = len(df_zmax)
global_pval = num_exceeding / num_toys
global_sig = stats.norm.ppf(1 - global_pval)

print(f"Global p-value (asymptotics): {100 * global_pval:.2f} %")
print(f"Global significance (asymptotics): {global_sig:.2f}")


# Calculate confidence intervals on global p and sig
global_pval_lo, global_pval_hi = binom_mle_interval(num_exceeding, num_toys)
print("68% CL interval on p-global: "
      f"[{100 * global_pval_lo:.2f} %, {100 * global_pval_hi:.2f} %]")

global_sig_lo, global_sig_hi = \
    stats.norm.ppf(1 - global_pval_hi), stats.norm.ppf(1 - global_pval_lo)
print("68% CL interval on sig-global: "
      f"[{global_sig_lo:.2f}, {global_sig_hi:.2f}]")


# Fit the trial factor
trial_factor = fit_trial_factor(df_zmax["max_sig"])
tf_global_sig = stats.norm.ppf(1 - (1 - (1 - obs_pval)**trial_factor))
print(f"Trial factor: {trial_factor:.1f}")
print(f"Global significance (trial factor): {tf_global_sig:.2f}")


# PDF for trial factor for plot
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
g_tf.Draw("C,SAME")

# Labels
latex.DrawLatex(0.54, 0.85, f"Number of toys: {len(df_zmax)}")
latex.DrawLatex(0.54, 0.80,
                "Toys with Z_{0}^{max} > Z_{0}^{max,obs}: "
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
latex.DrawLatex(0.54, 0.56, "Trial factor: "
                f"{trial_factor:.1f} (Z_{{global}} = {tf_global_sig:.2f})")

# Line at observed significance
line = R.TLine(obs_z0, 0, obs_z0, h_zmax.GetMaximum())
line.SetLineStyle(R.kDashed)
line.Draw()
latex.SetNDC(0)

latex.DrawLatex(obs_z0 + 0.05, 0.2 * h_zmax.GetMaximum(),
                "Z_{{0}}^{{max,obs}}(asymptotics) = {:.2f}".format(obs_z0))

c.RedrawAxis()
c.SaveAs("zmax_asymptotics.pdf")
