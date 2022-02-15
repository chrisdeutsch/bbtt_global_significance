#!/usr/bin/env python
import argparse
import sys
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("-o", "--outfile", default=None)
args = parser.parse_args()


# Read global significance toys
df = pd.read_csv(args.infile)

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

# Failed fits
df["failed_fit"] = (df["uncond_status"] != 0) | (df["cond_status"] != 0)

# Print or write failed fits to file for re-running
df_failed = df.loc[df["failed_fit"], ["toyindex", "mass"]]

if args.outfile is not None:
    df_failed.to_csv(args.outfile, index=False, header=False)
else:
    df_failed.to_csv(sys.stdout, index=False, header=False)
