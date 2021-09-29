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

parser.add_argument("--mu-range", default=15., type=float)
parser.add_argument("--optimizer-strategy", type=int, default=2)

parser.add_argument("-o", "--outfile", default=None)
parser.add_argument("-m", "--mass", type=int, default=None)
parser.add_argument("-i", "--index", type=int, default=None)

parser.add_argument("-v", "--verbose", action="store_true")

parser.add_argument("--globs-tree", default="",
                    help="Tree containing the global observables")
parser.add_argument("--globs-index", default=0, type=int,
                    help="Index in the tree that contains the values of the global observables")

args = parser.parse_args()

# Get macro directory to load needed ROOT macros
macro_path = os.path.join(os.path.dirname(__file__), "..", "macros")
macro_path = os.path.abspath(macro_path)


import ROOT as R
R.gROOT.SetBatch(True)
R.gROOT.ProcessLine(".L {}/DiscoveryTestStat.C+".format(macro_path))

R.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")
R.Math.MinimizerOptions.SetDefaultStrategy(args.optimizer_strategy)

ret = R.DiscoveryTestStat(
    args.infile,
    args.workspace_name,
    args.model_config,
    args.data_name,
    args.mu_range,
    args.globs_tree,
    args.globs_index,
    args.verbose)

# Warning: test statistic is the likelihood ratio and not q0: q0 = 2 * LLR
print("Likelihood-ratio: {:.5f}".format(ret.ts))
print("q0: {:.5f}".format(2 * ret.ts))
print("muhat: {:.5f}".format(ret.muhat))
print("uncond_status: {}".format(ret.uncond_status))
print("uncond_minNLL: {}".format(ret.uncond_minNLL))
print("cond_status: {}".format(ret.cond_status))
print("cond_minNLL: {}".format(ret.cond_minNLL))
print("cond_zhf: {}".format(ret.cond_zhf))
print("uncond_zhf: {}".format(ret.uncond_zhf))
print("cond_ttbar: {}".format(ret.cond_ttbar))
print("uncond_ttbar: {}".format(ret.uncond_ttbar))


if args.outfile:
    with open(args.outfile, "w") as f:
        fieldnames = [
            "index", "mass", "q0", "muhat",
            "uncond_status", "uncond_minNLL",
            "cond_status", "cond_minNLL",
            "cond_zhf", "uncond_zhf",
            "cond_ttbar", "uncond_ttbar",
            "mu_range"
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
            "cond_zhf": ret.cond_zhf,
            "uncond_zhf": ret.uncond_zhf,
            "cond_ttbar": ret.cond_ttbar,
            "uncond_ttbar": ret.uncond_ttbar,
            "mu_range": args.mu_range,
        })
