#!/usr/bin/env python
import argparse
import os
import re

import numpy as np
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("infiles", nargs="+")
parser.add_argument("-o", "--outfile", default=None)
args = parser.parse_args()


expected_seeds = set(range(1000))
pattern = re.compile(r"^toys_combined_(\d+)\.csv$")

retry_jobs = []

for fn in args.infiles:
    basename = os.path.basename(fn)

    m = pattern.match(basename)
    if not m:
        raise RuntimeError("Cannot parse filename")

    mass, = m.groups()
    mass = int(mass)

    df = pd.read_csv(fn)

    seeds = df["seed"].unique()
    finished_jobs = len(seeds)
    missing_seeds = expected_seeds - set(seeds)
    num_missing = len(missing_seeds)

    retry_jobs += [(mass, seed) for seed in missing_seeds]

    print(f"For file: {basename}")
    print(f"Number of finished jobs: {finished_jobs}")
    print(f"Number of missing jobs: {num_missing}")
    print(f"")

# Write out csv of jobs to rerun
if args.outfile is not None:
    num_retry = len(retry_jobs)
    print(f"Number of jobs to retry: {num_retry}")
    print(f"Writing to: {args.outfile}")

    with open(args.outfile, "w") as fout:
        for mass, seed in retry_jobs:
            fout.write(f"{mass},{seed}\n")
