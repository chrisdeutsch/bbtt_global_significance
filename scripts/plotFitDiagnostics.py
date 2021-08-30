#!/usr/bin/env python
import argparse
import re
import numpy as np
import pandas as pd

from os import path

parser = argparse.ArgumentParser()
parser.add_argument("infiles", nargs="+")
args = parser.parse_args()

import ROOT as R
R.gROOT.SetBatch(True)
R.gROOT.SetStyle("ATLAS")


dfs = []

# Guess masses from filename
for infile in args.infiles:
    basename = path.basename(infile)

    match = re.search(r"(\d{3,})", basename)
    assert match is not None

    mass = int(match.group(0))

    print("Got file '{}' with mass {}".format(basename, mass))

    df = pd.read_csv(infile)
    df["masspoint"] = mass
    dfs.append(df)

df = pd.concat(dfs)
df.status_uncond = df.status_uncond.astype(np.int64)
df.status_cond = df.status_cond.astype(np.int64)
df.rename(columns={"index": "toyindex"}, inplace=True)

df["fit_failed"] = (df.status_uncond != 0) | (df.status_cond != 0)
df["fit_success"] = ~df.fit_failed

# Average time per fit
avg_time = df.groupby("masspoint")["avg_time"].mean()
# Fraction of failed fits
failure_rate = df.groupby("masspoint")["fit_failed"].mean()
#
max_muhat = df.groupby("masspoint")["muhat"].max()
print(max_muhat)

p_success_all = df.groupby("masspoint")["fit_success"].mean().product()
print("Probability of having a successfully fitted toy for all points: {:.1f} %".format(100 * p_success_all))

# Average time histogram / barchart
h_avg_time = R.TH1F("h_avg_time", "", len(avg_time), 0, len(avg_time))
for idx, (mass, time) in enumerate(zip(avg_time.index, avg_time), start=1):
    h_avg_time.SetBinContent(idx, time)
    h_avg_time.GetXaxis().SetBinLabel(idx, "{}".format(mass))

h_avg_time.SetMinimum(0)
h_avg_time.SetMaximum(1.2 * avg_time.max())
h_avg_time.GetXaxis().SetTitle("Masspoint [GeV]")
h_avg_time.GetYaxis().SetTitle("Time per fit [s]")
h_avg_time.GetXaxis().SetLabelFont(43)
h_avg_time.GetXaxis().SetLabelSize(14)

# Failure rate histogram / barchart
h_failure_rate = R.TH1F("h_failure_rate", "", len(failure_rate), 0, len(failure_rate))
for idx, (mass, rate) in enumerate(zip(failure_rate.index, failure_rate), start=1):
    h_failure_rate.SetBinContent(idx, 100 * rate)
    h_failure_rate.GetXaxis().SetBinLabel(idx, "{}".format(mass))

h_failure_rate.SetMinimum(0)
h_failure_rate.SetMaximum(1.2 * 100 * failure_rate.max())
h_failure_rate.GetXaxis().SetTitle("Masspoint [GeV]")
h_failure_rate.GetYaxis().SetTitle("Fit Failure Rate [%]")
h_failure_rate.GetXaxis().SetLabelFont(43)
h_failure_rate.GetXaxis().SetLabelSize(14)

c = R.TCanvas("c", "", 800, 600)
h_avg_time.Draw("HIST")
c.SaveAs("fit_avg_time.pdf")

c.Clear()
h_failure_rate.Draw("HIST")
c.SaveAs("fit_failure_rate.pdf")
