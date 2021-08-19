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
