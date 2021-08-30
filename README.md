# Estimation of Global Significance for HH->bbtautau

ROOT version used:
`lsetup "root 6.20.06-x86_64-centos7-gcc8-opt"`

## Brief Introduction

`runDiscoveryTestStat.py`:

Runs the fits to determine the discovery test statistic q0 on
data. Usually the data is replaced with pseudodata generated in
CxAODReader. This is used to estimate the effect of the look-elsewhere
effect.


`runDiscoveryTestStatToys.py`:

Generates toys from a given workspace under the b-only hypothesis and
determines the discovery test statistic q0. It stores the results for
every toy in an output `.csv` file. Afterwards multiple runs of this
tool (with different seeds) can be combined to get the q0 sampling
distribution under the b-only hypothesis.


`plotFitDiagnostics.py`:

Creates a couple of diagnostic plots (avg. time per toy, fit failure
rate) from a set of toys.


## Batch Submission

The scripts in `batch_submission` are written to be run from inside of
this directory and 'hardcoded' to work in Bonn. Please adapt the
scripts before using it elsewhere. They should however serve as a
reasonable starting point.

Before submitting make sure the relevant directories exist:
```bash
mkdir logs
mkdir /cephfs/user/s6crdeut/WSMakerToys_q0/combined_{251,260,280,300,325,350,375,400,450,500,550,600,700,800,900,1000,1100,1200,1400,1600}
```

To submit the q0 sampling distribution toys for the `mX = 500 GeV` workspace:
```bash
condor_submit Mass=500 SeedOffset=0 submission_toys_local.jdl
```

There is a convenient script `run_local_toys.sh` that runs all mass
points (on batch) successively and merges the results.


## Evaluation of Results
