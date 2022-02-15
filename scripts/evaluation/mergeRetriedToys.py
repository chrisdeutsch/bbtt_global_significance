#!/usr/bin/env python
import argparse

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("outfile")
parser.add_argument("toys")
parser.add_argument("retries")
args = parser.parse_args()


# Toys before retrying
df = pd.read_csv(args.toys)

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

# Check that there are no duplicates
assert not df.duplicated(["toyindex", "mass"]).any()

# Flag failed fits
df["failed_fit"] = (df["uncond_status"] != 0) | (df["cond_status"] != 0)

# Load retries and merge with original toys
df_retried = pd.read_csv(args.retries)

# Check that there are no duplicates
assert not df_retried.duplicated(["toyindex", "mass"]).any()

# Merge original toys with retries dropping failed fits
df_merged = pd.concat([df, df_retried])
df_merged.sort_values("failed_fit", inplace=True)
df_merged.drop_duplicates(["toyindex", "mass"], keep="first", inplace=True)

# Ensure that all toys have 20 fits
assert np.all(df_merged.groupby("toyindex")["q0"].count() == 20)

# Store that puppy
df_merged.to_csv(args.outfile, index=False)
