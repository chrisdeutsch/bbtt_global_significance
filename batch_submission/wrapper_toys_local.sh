#!/usr/bin/env bash
(($# == 4)) || { echo "wrapper_toys_local.sh infile outdir ntoys seed"; exit 1; }

infile="${1}"
outdir="${2}"
ntoys="${3}"
seed="${4}"

[[ -f "${infile}" ]] || { echo "Infile ${infile} does not exist"; exit 1; }
[[ -d "${outdir}" ]]  || { echo "Outdir does not exist"; exit 1; }

mkdir -p /jwd/run
pushd /jwd/run

# Get code
tar -xzf /cephfs/user/s6crdeut/bbtt_global_significance.tar.gz \
    || { echo "Cannot get code"; exit 1; }

if [ -d /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase ]; then
    export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
    source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
else
    exit 1
fi

# Recommended CentOS7 root version
lsetup "root 6.20.06-x86_64-centos7-gcc8-opt"

# Set up PATH
script_dir="$(readlink -e bbtt_global_significance/scripts)"
PATH="${script_dir}:${PATH}"

outfile="${outdir}/toy_${seed}.csv"
runDiscoveryTestStatToys.py "${infile}" -s "${seed}" -n "${ntoys}" -o "${outfile}" 2>&1
