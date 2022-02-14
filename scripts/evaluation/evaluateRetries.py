#!/usr/bin/env python
from glob import glob
import argparse
import numpy as np
import os
import pandas as pd
import re
import tarfile


parser = argparse.ArgumentParser()
parser.add_argument("indir")
parser.add_argument("-o", "--outfile", default=None)
args = parser.parse_args()


pattern = re.compile(r"^retry_(.*)_idx\d+_m\d+.csv$")


def get_retry_df(fn):
    dfs = []
    with tarfile.open(fn, "r") as tar:
        for name in tar.getnames():
            m = pattern.match(name)
            if not m:
                raise RuntimeError("Cannot parse: " + name)

            retry_method, = m.groups()

            df = pd.read_csv(tar.extractfile(name))
            df["retry_method"] = retry_method
            dfs.append(df)

    return pd.concat(dfs)


# Helper to print results nicely
def print_toy(df, index, mass):
    print_cols = ["toyindex", "mass",
                  "q0", "muhat", "mu_range",
                  "uncond_minNLL", "cond_minNLL",
                  "uncond_status", "cond_status",
                  "retry_method", "retry_priority"]

    print(df.loc[
        (df["toyindex"] == index) & (df["mass"] == mass), print_cols]
          .sort_values("retry_priority"))


dfs = []
for fn in glob(os.path.join(args.indir, "retry_*.tar.gz")):
    df = get_retry_df(fn)
    dfs.append(df)

# Number of retried toys
num_retried = len(dfs)
print(f"Number of retried toys: {num_retried}")

# Make joint dataframe
df = pd.concat(dfs)

# Highest (0) to lowest (6)
priority = {
    "minuit2strat2": 0,
    "minuit2strat1mu2": 1,
    "minuit2strat1mu0p5": 2,
    "minuitstrat2": 3,
    "minuitstrat1": 4,
    "minuitstrat1mu2": 5,
    "minuitstrat1mu0p5": 6,
}

df["retry_priority"] = df["retry_method"].map(priority)


# Convert to proper dtypes
df = df.astype({
    "uncond_status": "int64",
    "cond_status": "int64",
    "uncond_covQual": "int64",
    "cond_covQual": "int64",
})

# Rename 'index' to 'toyindex' to avoid confusion with the index of
# the dataframe
df.rename(columns={"index": "toyindex"}, inplace=True)

# Filter out bad fits
df["failed_fit"] = (df["uncond_status"] != 0) | (df["cond_status"] != 0)
df["neg_q0"] = (df["q0"] < -1e-1)

sel_good = (~df["failed_fit"]) & (~df["neg_q0"])
df_good = df.loc[sel_good].copy()

# Count number of successful fits
df_good["good_count"] = df_good.groupby(["toyindex", "mass"])["failed_fit"] \
                               .transform("count")

# Require at least three good alternative fits
df_good = df_good.loc[df_good["good_count"] >= 3].copy()


# === Sanity check ===
# Check the spread of q0/sig for good fits
df_good_copy = df_good.copy()
df_good_copy["sig"] = np.sqrt(np.abs(df_good_copy["q0"]))

df_sanity = df_good_copy.groupby(["toyindex", "mass"]).agg(
    q0_max=("q0", "max"),
    q0_med=("q0", "median"),
    q0_min=("q0", "min"),
    sig_max=("sig", "max"),
    sig_med=("sig", "median"),
    sig_min=("sig", "min"),
    count=("q0", "count")
)

df_sanity["sig_delta"] = df_sanity["sig_max"] - df_sanity["sig_min"]

print("Suspicious fits (included in result):")
with pd.option_context("display.max_rows", None,
                       "display.max_columns", None,
                       "display.width", 185,
                       "display.precision", 2):
    print(df_sanity.loc[df_sanity["sig_delta"] > 0.2]
          .sort_values("sig_delta", ascending=False))

# === End of sanity check ===

# Sort by priority
df_good.sort_values("retry_priority", inplace=True)

# Remove duplicates (keep first -> highest priority)
df_good.drop_duplicates(["toyindex", "mass"], keep="first", inplace=True)

# Good fits after retry
num_good = len(df_good)
print(f"Good toys after retrying: {num_good}")

print("Count of chosen retry method:")
print(df_good.groupby("retry_method")["toyindex"].count())

if args.outfile is not None:
    df_good.to_csv(args.outfile, index=False)


df_bad = df.merge(df_good,
                  on=["toyindex", "mass"],
                  how="outer",
                  suffixes=("", "_other"),
                  indicator=True)
df_bad = df_bad.loc[df_bad["_merge"] == "left_only"]

df_bad["uncond_good"] = df_bad["uncond_status"] == 0
df_bad["cond_good"] = df_bad["cond_status"] == 0
df_bad["muhat_negative"] = df_bad["muhat"] < 0

print("Count of status codes")
print(df_bad.groupby(["uncond_good", "cond_good"])
            .agg(count=("q0", "count"), muhat_sign=("muhat_negative", "mean")))

print(df_bad.groupby(["uncond_good", "cond_good", "muhat_negative"])["q0"]
            .count())
