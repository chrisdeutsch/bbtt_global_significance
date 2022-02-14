#!/usr/bin/env python
import argparse
import re
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

from scipy import stats
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import load_toys, ToyPvalueCalculator


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
          "-- Setting significance to 4.99")
    df_good.loc[mask_inf, "sig"] = 4.99


# Calculate global significance
df_zmax = df_good.groupby("toyindex")["sig"].max().to_frame("max_sig")

num_exceeding = (df_zmax["max_sig"] > obs_z0).sum()
num_toys = len(df_zmax)
global_pval = num_exceeding / num_toys
global_sig = stats.norm.ppf(1 - global_pval)

print(f"Global p-value: {global_pval:.4f}")
print(f"Global significance: {global_sig:.4f}")


# Bootstrapping
num_bootstraps = 1000
zglobal_bootstraps = []
rng = np.random.default_rng(1234597654)

for i in tqdm(range(num_bootstraps), "Bootstrapping"):
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
        df_good.loc[mask_inf, "sig"] = 4.99

    # Calculate global significance
    df_zmax = df_good.groupby("toyindex")["sig"].max().to_frame("max_sig")

    num_exceeding = (df_zmax["max_sig"] > obs_z0).sum()
    num_toys = len(df_zmax)
    global_pval = num_exceeding / num_toys
    global_sig = stats.norm.ppf(1 - global_pval)

    zglobal_bootstraps.append(global_sig)

zglobal_bootstraps = np.array(zglobal_bootstraps)

print("Global significance (bootstrap mean): "
      f"{zglobal_bootstraps.mean():.3f}")
print("Global significance error (bootstrap std): "
      f"{zglobal_bootstraps.std():.3f}")
print("Error on mean global significance: "
      f"{zglobal_bootstraps.std() / np.sqrt(len(zglobal_bootstraps)):.3f}")


fig, ax = plt.subplots()
ax.hist(zglobal_bootstraps, range=(1.8, 2.1), bins=50)
ax.set_xlim(r"$Z_\mathrm{global}$")
ax.set_ylim(r"Bootstrap Samples")
fig.savefig("test.pdf")
