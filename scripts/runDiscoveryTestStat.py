#!/usr/bin/env python
import argparse
import csv
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("--workspace-name", default="combined")
parser.add_argument("--model-config", default="ModelConfig")
parser.add_argument("--data-name", default="obsData")

parser.add_argument("-o", "--outfile", default=None)
parser.add_argument("-m", "--mass", type=int, default=None)
parser.add_argument("-i", "--index", type=int, default=None)

args = parser.parse_args()

# Get macro directory to load needed ROOT macros
macro_path = os.path.join(os.path.dirname(__file__), "..", "macros")
macro_path = os.path.abspath(macro_path)


import ROOT as R
R.gROOT.SetBatch(True)
R.gROOT.ProcessLine(".L {}/DiscoveryTestStat.C++".format(macro_path))

R.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")
R.Math.MinimizerOptions.SetDefaultStrategy(2)

ret = R.DiscoveryTestStat(
    args.infile,
    args.workspace_name,
    args.model_config,
    args.data_name)


# Warning: test statistic is the likelihood ratio and not q0: q0 = 2 * LLR
print("Likelihood-ratio: {:.5f}".format(ret.ts))
print("q0: {:.5f}".format(2 * ret.ts))
print("muhat: {:.5f}".format(ret.muhat))
print("uncond_status: {}".format(ret.uncond_status))
print("uncond_minNLL: {}".format(ret.uncond_minNLL))
print("cond_status: {}".format(ret.cond_status))
print("cond_minNLL: {}".format(ret.cond_minNLL))


if args.outfile:
    with open(args.outfile, "w") as f:
        fieldnames = [
            "index", "mass", "q0", "muhat",
            "uncond_status", "uncond_minNLL",
            "cond_status", "cond_minNLL"
        ]

        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        writer.writerow({
            "index": args.index,
            "mass": args.mass,
            "q0": 2 * ret.ts,
            "muhat": ret.muhat,
            "uncond_status": ret.uncond_status,
            "uncond_minNLL": ret.uncond_minNLL,
            "cond_status": ret.cond_status,
            "cond_minNLL": ret.cond_minNLL,
        })
