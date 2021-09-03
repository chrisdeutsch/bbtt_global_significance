#!/usr/bin/env python
import argparse
import re
import os
from os import path
import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.optimize import root_scalar

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("--method", choices=["asymptotics", "toys"], default="asymptotics")
parser.add_argument("-o", "--outdir", default="")
args = parser.parse_args()


import ROOT as R
R.gROOT.SetBatch(True)
R.gROOT.SetStyle("ATLAS")


# CONFIGURE ME!
obs_z0 = 3.04846308291
obs_pval = 1.0 - stats.norm.cdf(obs_z0)


# Read inputs and perform preprocessing
df = pd.read_csv(args.infile)
df["uncond_status"] = df["uncond_status"].astype(np.int)
df["cond_status"] = df["cond_status"].astype(np.int)
df["mass"] = df["mass"].astype(np.int)
df.rename(columns={"index": "toyindex"}, inplace=True)
print("Total number of toys: {}".format(len(np.unique(df["toyindex"]))))


# Check for failed fits
df["failed_fit"] = (df["uncond_status"] != 0) | (df["cond_status"] != 0)
num_failed = df["failed_fit"].sum()
frac_failed = num_failed / float(len(df))
print("Failed fits: {} ({:.1f} %)".format(num_failed, 100 * frac_failed))

failed_fits = df.loc[df["failed_fit"]].copy()
failed_fits["muhatPositive"] = failed_fits["muhat"] >= 0.0
print(failed_fits.groupby("muhatPositive").count())

# # Check what happens when we accept fits with IMPROVE errors when muhat is negative
# uncond_good = (df["uncond_status"] == 0) | (df["uncond_status"] == 4000)
# cond_good = (df["cond_status"] == 0) | (df["cond_status"] == 4000)
# neg_muhat = df["muhat"] < 0
# df.loc[df["failed_fit"] & cond_good & uncond_good & neg_muhat, "q0"] = 0
# df.loc[cond_good & uncond_good & neg_muhat, "failed_fit"] = False

# =====================
# Basic fit diagnostics
# =====================
fit_df = df.loc[~df["failed_fit"]] \
           .groupby("mass")["muhat"] \
           .agg([("muhatMin", "min"),
                 ("muhatMax", "max"),
                 ("muhatMean", "mean"),
                 ("muhatStd", "std"),
                 ("nToys", "count")])
fit_df["muhatMeanError"] = fit_df["muhatStd"] / np.sqrt(fit_df["nToys"])
fit_df["muBiasSig"] = fit_df["muhatMean"] / fit_df["muhatMeanError"]

mu_range = df.loc[~df["failed_fit"]] \
             .groupby("mass")["mu_range"].min().rename("muRange")
fit_df = pd.concat([fit_df, mu_range], axis=1)
fit_df["maxMuhatOverRange"] = fit_df["muhatMax"] / fit_df["muRange"]

rate_fit_fail = df.groupby("mass")["failed_fit"].mean().rename("rateFitFail")
fit_df = pd.concat([fit_df, rate_fit_fail], axis=1)

fit_df["muhatAbsMax"] = np.maximum(-fit_df["muhatMin"], fit_df["muhatMax"])
fit_df["muhatAbsMaxOverStd"] = fit_df["muhatAbsMax"] / fit_df["muhatStd"]

print(fit_df)

# Muhat plots to check for bias
fig, axs = plt.subplots(nrows=5, ncols=4, sharey=True, figsize=(14, 12))
for ax, (mass, mass_df) in zip(axs.ravel(), df.loc[~df["failed_fit"]].groupby("mass")):
    std = mass_df["muhat"].std()
    mean = mass_df["muhat"].mean()

    ax.hist(mass_df["muhat"], histtype="step", bins=25, range=(-5 * std, 5 * std))
    ax.set_xlabel(r"$\hat{\mu}$")
    ax.set_ylabel(r"Toy Experiments")

    ax.annotate("$m_{{X}}$ = {} GeV".format(mass), xy=(0.57, 0.85), xycoords="axes fraction")
    ax.annotate(r"$\langle \hat{\mu} \rangle$ = "
                + "{:.2f}".format(mean) # Mean
                + r"$\pm$"
                + "{:.2f}".format(std / np.sqrt(len(mass_df))), # Error of mean
                xy=(0.57, 0.7), xycoords="axes fraction")

fig.tight_layout()
fig.savefig(path.join(args.outdir, "muhat.pdf"))

# Need to require that all fits were successful for a given toy
df["good_toy"] = ~df.groupby("toyindex")["failed_fit"].transform(np.any)
print("Fraction of good toys: {:.1f} %".format(100 * df["good_toy"].sum() / float(len(df))))

# Fix q0 for muhat <= 0 for one-sided discovery test statistic
df.loc[df["muhat"] <= 0, "q0"] = 0.0

# Select out only good toys
good_df = df.loc[df["good_toy"]].copy()

# Set slightly negative q0 to 0
threshold = -0.05
good_df.loc[(good_df["q0"] < 0) & (good_df["q0"] > threshold), "q0"] = 0.0

# Hopefully there are no more negative q0's left
#print(good_df.loc[good_df["q0"] < 0])
assert len(good_df.loc[good_df["q0"] < 0]) == 0

# Significance using asymptotic approximation
good_df["sig_asymptotics"] = np.sqrt(good_df["q0"])
good_df["sig_toys"] = 0.0 # TODO

zmax_df = good_df.groupby("toyindex")["sig_asymptotics"].max().to_frame("max_sig_asymptotics")
zmax_df["max_sig_toys"] = good_df.groupby("toyindex")["sig_toys"].max()


# Calculate global p-value / significance
num_exceeding = (zmax_df["max_sig_asymptotics"] > obs_z0).sum()
num_toys = len(zmax_df)
global_pval = num_exceeding / float(num_toys)
global_sig = stats.norm.ppf(1.0 - global_pval)


# Uncertainties
llh = lambda p, k, n: k * R.TMath.Log(p) + (n - k) * R.TMath.Log(1 - p)
delta_llh = lambda p: llh(p, num_exceeding, num_toys) - llh(global_pval, num_exceeding, num_toys)

root_lo = root_scalar(lambda x: -2 * delta_llh(x) - 1, bracket=[0.01, global_pval])
root_hi = root_scalar(lambda x: -2 * delta_llh(x) - 1, bracket=[global_pval, 0.1])
assert root_lo.converged and root_hi.converged

global_pval_up = root_hi.root - global_pval
global_pval_dn = global_pval - root_lo.root

global_sig_up = stats.norm.ppf(1.0 - root_lo.root) - global_sig
global_sig_dn = global_sig - stats.norm.ppf(1.0 - root_hi.root)


# ==========
#  Plotting
# ==========

# Maximum local significance plot
h_zmax = R.TH1F("h_zmax", "", 120, 0, 6)
for index, row in zmax_df.iterrows():
    h_zmax.Fill(row["max_sig_asymptotics"])

h_zmax.GetXaxis().SetTitle("Maximum local significance Z_{0}^{max}")
h_zmax.GetYaxis().SetTitle("Toy experiments")
h_zmax.SetMinimum(0)
h_zmax.SetMaximum(1.3 * h_zmax.GetMaximum())
h_zmax.SetLineColor(R.kBlue)
h_zmax.SetMarkerSize(0)

latex = R.TLatex()
latex.SetNDC()
latex.SetTextFont(43)
latex.SetTextSize(19)

c = R.TCanvas("c", "", 800, 600)
h_zmax.Draw("HIST,E0")

# Labels
latex.DrawLatex(
    0.54, 0.85,
    "Number of toys: {}".format(len(zmax_df)))

latex.DrawLatex(
    0.54, 0.80,
    "Toys with Z_{0}^{max} > Z_{0}^{max,obs}: " + "{}".format(num_exceeding))

latex.DrawLatex(
    0.54, 0.72,
    "Global p-value (toys): {:.4f}#splitline{{#plus {:.4f}}}{{#minus {:.4f}}}".format(
        global_pval, global_pval_up, global_pval_dn))

latex.DrawLatex(
    0.54, 0.64,
    "Global significance (toys): {:.2f}#splitline{{#plus {:.2f}}}{{#minus {:.2f}}}".format(
        global_sig, global_sig_up, global_sig_dn))

# Observed max. significance line
line = R.TLine(obs_z0, 0, obs_z0, h_zmax.GetMaximum())
line.SetLineStyle(R.kDashed)
line.Draw()

latex.SetNDC(0)
latex.DrawLatex(2.83 + 0.05, 0.4 * h_zmax.GetMaximum(),
                    "Z_{{0}}^{{max,obs}}({}) = {:.2f}".format(args.method, obs_z0))

c.RedrawAxis()
c.SaveAs(path.join(args.outdir, "zmax.pdf"))


# Fit failure rate plot
h_failrate = R.TH1F("h_failrate", "", len(fit_df), 0, len(fit_df))
for idx, (mass, rate) in enumerate(fit_df["rateFitFail"].items(), start=1):
    h_failrate.SetBinContent(idx, 100 * rate)
    h_failrate.GetXaxis().SetBinLabel(idx, "{}".format(mass))

h_failrate.SetMinimum(0)
h_failrate.SetMaximum(1.2 * 100 * fit_df["rateFitFail"].max())
h_failrate.GetXaxis().SetTitle("Masspoint [GeV]")
h_failrate.GetYaxis().SetTitle("Rate of Failing Fits [%]")
h_failrate.GetXaxis().SetLabelFont(43)
h_failrate.GetXaxis().SetLabelSize(14)

c.Clear()
h_failrate.Draw("HIST")
c.RedrawAxis()
c.SaveAs(path.join(args.outdir, "failrate.pdf"))
