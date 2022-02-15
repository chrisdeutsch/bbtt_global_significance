#!/usr/bin/env python
import argparse
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

from scipy.stats import norm
from statsmodels.distributions.empirical_distribution import ECDF
from tqdm import tqdm
import numpy as np
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("-m", "--mass", default="251")
parser.add_argument("--outfile-q0-plot", default="q0.pdf")
parser.add_argument("--outfile-sig-plot", default="sig.pdf")
parser.add_argument("--outfile-q0-csv", default="q0.csv")
args = parser.parse_args()


import ROOT as R
R.gROOT.SetBatch(True)
R.gROOT.SetStyle("ATLAS")


# Read local significance toys
df = pd.read_csv(args.infile)

# Convert to proper dtypes
df = df.astype({
    "uncond_status": "int64",
    "cond_status": "int64",
    "uncond_covQual": "int64",
    "cond_covQual": "int64",
})

# Check that there are no duplicates
assert not df.duplicated(["seed", "index"]).any()

# Failed fits
df["failed_fit"] = (df["uncond_status"] != 0) | (df["cond_status"] != 0)
df.loc[df["failed_fit"], "q0"] = 0.0

num_failed = df["failed_fit"].sum()
frac_failed = num_failed / len(df)
print(f"Failed fits: {num_failed} ({100 * frac_failed:.1f} %)")


# Transform q0 to one-sided discovery test statistic
df.loc[df["muhat"] <= 0, "q0"] = 0.0

# Only keep good toys
df_good = df.loc[~df["failed_fit"]].copy()

# Set negative q0 to 0
df_good.loc[df_good["q0"] < 0, "q0"] = 0.0


# Asymptotic approximation
# 1/2 delta(q0) + 1/2 chi^2(q0, NDF=1) density histogram
h_q0_sampling = R.TH1F("h_q0_sampling", "", 40, 0, 20)
h_q0_sampling.SetLineColor(R.kRed)
h_q0_sampling.SetMarkerColor(R.kRed)
for bin_idx in range(1, h_q0_sampling.GetNbinsX() + 1):
    xaxis = h_q0_sampling.GetXaxis()

    hi_edge = xaxis.GetBinUpEdge(bin_idx)
    lo_edge = xaxis.GetBinLowEdge(bin_idx)
    width = hi_edge - lo_edge

    # Chi-squared part of the density
    proba = 0.5 * (R.Math.chisquared_cdf(hi_edge, 1, 0)
                   - R.Math.chisquared_cdf(lo_edge, 1, 0))

    # Delta function part of the density
    eps = 1e-9
    if lo_edge <= eps and eps < hi_edge:
        proba += 0.5

    h_q0_sampling.SetBinContent(bin_idx, proba / width)


# q0 toy histogram
h_q0 = R.TH1F("h_q0", "", 40, 0, 20)
h_q0.FillN(len(df_good), df_good["q0"].values, R.nullptr)
h_q0.GetXaxis().SetTitle("q_{0}")
h_q0.GetYaxis().SetTitle("Probability Density")
h_q0.Scale(1.0 / len(df_good), "WIDTH")
h_q0.SetLineColor(R.kBlue)
h_q0.SetMarkerSize(0)
h_q0.SetMinimum(1e-5)


leg = R.TLegend(0.6, 0.72, 0.85, 0.85)
leg.SetTextFont(43)
leg.SetTextSize(19)
leg.SetBorderSize(0)
leg.AddEntry(h_q0, "Toys", "f")
leg.AddEntry(h_q0_sampling, "Asymptotic approx.", "f")

latex = R.TLatex()
latex.SetTextFont(43)
latex.SetTextSize(19)
latex.SetNDC()

c = R.TCanvas("c", "", 800, 600)
c.SetLogy()

h_q0.Draw("HIST,E0")
h_q0_sampling.Draw("HIST,SAME")
leg.Draw()

latex.DrawLatex(0.55, 0.65,
                "Combined #tau_{had}#tau_{had}, "
                "#tau_{lep}#tau_{had} (SLT, LTT)")
latex.DrawLatex(0.55, 0.6, f"m_{{X}} = {args.mass} GeV")

c.RedrawAxis()
c.SaveAs(args.outfile_q0_plot)
del c


# ECDF
cdf = ECDF(df_good["q0"])

# Store q0 values for global significance analysis
df_good["q0"].to_csv(args.outfile_q0_csv, index=False)

# Significance plot
q0 = np.linspace(0, 14, 200)
sig = norm.ppf(cdf(q0))

# Bootstrap uncertainty on sig
band = []
for i_bootstrap in tqdm(range(1000), "Bootstrap"):
    bootstrap_sample = df_good["q0"].sample(frac=1, replace=True)
    bootstrap_cdf = ECDF(bootstrap_sample)
    band.append(norm.ppf(bootstrap_cdf(q0)))

p = (100 - 68) / 2
band = np.array(band)
band = np.percentile(band, [p, 100 - p], axis=0)

f = R.TF1("f", "sqrt(x)", 0, 14)
g = R.TGraph(len(q0), q0, sig)

err_hi = band[1] - sig
err_lo = sig - band[0]
g_err = R.TGraphAsymmErrors(len(q0), q0, sig,
                            R.nullptr, R.nullptr,
                            err_lo, err_hi)

f.SetLineColor(R.kRed)
f.SetLineWidth(2)
f.SetLineStyle(R.kDashed)

g.SetLineColor(R.kBlue)
g.SetLineWidth(2)

g_err.SetFillColor(R.kBlue - 9)

leg = R.TLegend(0.5, 0.4, 0.8, 0.55)
leg.SetTextFont(43)
leg.SetTextSize(19)
leg.SetBorderSize(0)
leg.AddEntry(g, "Toys", "l")
leg.AddEntry(g_err, "Uncertainty (68% CI)", "f")
leg.AddEntry(f, "Asymptotic approx.", "l")


h_dummy = R.TH1F("h_dummy", "", 5, 0, 14)
h_dummy.SetMinimum(0)
h_dummy.SetMaximum(4)
h_dummy.GetXaxis().SetTitle("q_{0}")
h_dummy.GetYaxis().SetTitle("Local Significance")

c = R.TCanvas("c", "", 800, 600)

h_dummy.Draw("AXIS")

g_err.Draw("E3,SAME")
g.Draw("L,SAME")
f.Draw("L,SAME")

leg.Draw()

latex.DrawLatex(0.5, 0.35,
                "Combined #tau_{had}#tau_{had}, "
                "#tau_{lep}#tau_{had} (SLT, LTT)")
latex.DrawLatex(0.5, 0.3, f"m_{{X}} = {args.mass} GeV")

c.RedrawAxis()
c.SaveAs(args.outfile_sig_plot)
