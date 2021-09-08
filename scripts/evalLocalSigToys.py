#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import re
from os import path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("infiles", nargs="+")
parser.add_argument("-o", "--outdir", default="")
args = parser.parse_args()


import ROOT as R
R.gROOT.SetBatch(True)
R.gROOT.SetStyle("ATLAS")


dfs = []
for infile in sorted(args.infiles):
    df = pd.read_csv(infile)

    # Guess masses from filename
    basename = path.basename(infile)

    match = re.search(r"(\d{3,})", basename)
    assert match is not None

    mass = int(match.group(0))

    print("Got file '{}' with mass {}".format(basename, mass))
    df["mass"] = mass

    dfs.append(df)

df = pd.concat(dfs)
df["status_uncond"] = df["status_uncond"].astype(np.int)
df["status_cond"] = df["status_cond"].astype(np.int)
df["mass"] = df["mass"].astype(np.int)
df.rename(columns={"index": "toyindex"}, inplace=True)

df["failed_fit"] = (df["status_uncond"] != 0) | (df["status_cond"] != 0)
num_failed = df["failed_fit"].sum()
frac_failed = num_failed / float(len(df))
print("Failed fits: {} ({:.1f} %)".format(num_failed, 100 * frac_failed))

# =====================
# Basic fit diagnostics
# =====================
fit_df = df.loc[~df["failed_fit"]] \
           .groupby("mass")["muhat"] \
           .agg([("muhatMean", "mean"),
                 ("muhatStd", "std"),
                 ("nToys", "count")])
fit_df["muhatMeanError"] = fit_df["muhatStd"] / np.sqrt(fit_df["nToys"])
fit_df["muBiasSig"] = fit_df["muhatMean"] / fit_df["muhatMeanError"]
fit_df["rateFitFail"] = df.groupby("mass")["failed_fit"].mean()
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
