#!/usr/bin/env python
import argparse
import csv
import os
import sys
import time

parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("-s", "--seed", type=int, required=True)
parser.add_argument("-o", "--outfile", default="toys.csv")
parser.add_argument("-n", "--ntoys", type=int, default=10)

parser.add_argument("--workspace-name", default="combined")
parser.add_argument("--model-config", default="ModelConfig")
parser.add_argument("--data-name", default="obsData")

parser.add_argument("--mu-range", default=15., type=float)

parser.add_argument("--optimizer-strategy", type=int, default=1)
parser.add_argument("--optimizer", choices=["Minuit2", "Minuit"], default="Minuit2")
parser.add_argument("-v", "--verbose", action="store_true")

args = parser.parse_args()

# Get macro directory to load needed ROOT macros
macro_path = os.path.join(os.path.dirname(__file__), "..", "macros")
macro_path = os.path.abspath(macro_path)

import ROOT as R
R.gROOT.SetBatch(True)
R.gROOT.ProcessLine(".L {}/DiscoveryTestStatToys.C++".format(macro_path))


# Retrieves value of a RooRealVar from RooArgSet
# Returns None if it does not exist
def retrieve_arg(argset, argname):
    if isinstance(argset.find(argname), R.RooRealVar):
        return argset.find(argname).getVal()

    return None


R.RooRandom.randomGenerator().SetSeed(10000 + args.seed)
R.Math.MinimizerOptions.SetDefaultMinimizer(args.optimizer)
R.Math.MinimizerOptions.SetDefaultStrategy(args.optimizer_strategy)

if args.verbose:
    # Doesn't really do anything...
    R.Math.MinimizerOptions.SetDefaultPrintLevel(3)

start_time = time.time()

htr = R.DiscoveryTestStatToys(
    args.infile,
    args.workspace_name,
    args.model_config,
    args.data_name,
    args.ntoys,
    args.mu_range,
    args.verbose)

end_time = time.time()

results = []
null_details = htr.GetNullDetailedOutput()
for i in range(null_details.numEntries()):
    argset = null_details.get(i)

    ts = retrieve_arg(argset, "ModelConfigB_only_TS0")
    muhat = retrieve_arg(argset, "ModelConfigB_only_TS0_fitUncond_SigXsecOverSM")

    uncond_status = retrieve_arg(argset, "ModelConfigB_only_TS0_fitUncond_fitStatus")
    uncond_minNLL = retrieve_arg(argset, "ModelConfigB_only_TS0_fitUncond_minNLL")

    cond_status = retrieve_arg(argset, "ModelConfigB_only_TS0_fitCond_fitStatus")
    cond_minNLL = retrieve_arg(argset, "ModelConfigB_only_TS0_fitCond_minNLL")

    cond_zhf = retrieve_arg(argset, "ModelConfigB_only_TS0_fitCond_ATLAS_norm_Zhf")
    uncond_zhf = retrieve_arg(argset, "ModelConfigB_only_TS0_fitUncond_ATLAS_norm_Zhf")

    cond_ttbar = retrieve_arg(argset, "ModelConfigB_only_TS0_fitCond_ATLAS_norm_ttbar")
    uncond_ttbar = retrieve_arg(argset, "ModelConfigB_only_TS0_fitUncond_ATLAS_norm_ttbar")

    uncond_covQual = retrieve_arg(argset, "ModelConfigB_only_TS0_fitUncond_covQual")
    cond_covQual = retrieve_arg(argset, "ModelConfigB_only_TS0_fitCond_covQual")

    results.append({
        "index": i,
        "seed": args.seed,
        "q0": 2 * ts,
        "muhat": muhat,
        "uncond_status": uncond_status,
        "uncond_minNLL": uncond_minNLL,
        "cond_status": cond_status,
        "cond_minNLL": cond_minNLL,
        "zhf_norm_cond": cond_zhf,
        "zhf_norm_uncond": uncond_zhf,
        "ttbar_norm_cond": cond_ttbar,
        "ttbar_norm_uncond": uncond_ttbar,
        "uncond_covQual": uncond_covQual,
        "cond_covQual": cond_covQual,
    })


total_time = end_time - start_time
time_per_toy = total_time / args.ntoys

print("Total time: {:2f} s".format(total_time))
print("Time per toy: {:2f} s/toy".format(time_per_toy))

with open(args.outfile, "w") as csvfile:
    fieldnames = [
        "q0", "muhat",
        "uncond_status", "uncond_minNLL",
        "cond_status", "cond_minNLL",
        "seed", "index",
        "avg_time", "mu_range",
        "zhf_norm_cond", "zhf_norm_uncond",
        "ttbar_norm_cond", "ttbar_norm_uncond",
        "uncond_covQual", "cond_covQual"]

    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in results:
        row["avg_time"] = time_per_toy
        row["mu_range"] = args.mu_range
        writer.writerow(row)
