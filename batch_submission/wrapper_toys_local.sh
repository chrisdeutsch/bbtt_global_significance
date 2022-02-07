#!/usr/bin/env bash
set -eu
(($# == 4)) || { echo "wrapper_toys_local.sh infile outdir ntoys seed"; exit 1; }

echo "Start: $(date)"

infile="${1}"
outdir="${2}"
ntoys="${3}"
seed="${4}"

[[ -f "${infile}" ]] || { echo "Infile ${infile} does not exist"; exit 1; }
[[ -d "${outdir}" ]]  || { echo "Outdir does not exist"; exit 1; }

mkdir -p /jwd/run
cd /jwd/run

# Get code
tar -xzf /cephfs/user/s6crdeut/bbtt_global_significance.tar.gz \
    || { echo "Cannot get code"; exit 1; }

set +eu
if [ -d /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase ]; then
    export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
    source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
else
    exit 1
fi

# Recommended CentOS7 root version
lsetup "root 6.20.06-x86_64-centos7-gcc8-opt"
set -eu

# Set range of mu
filename=$(basename "${infile}")
mass="${filename%.root}"

mu_range=""
case "${mass}" in
    251)  mu_range="10.0"  ;;
    260)  mu_range="18.0" ;;
    280)  mu_range="18.0" ;;
    300)  mu_range="14.0" ;;
    325)  mu_range="13.0" ;;
    350)  mu_range="8.0"  ;;
    375)  mu_range="5.0"  ;;
    400)  mu_range="4.0"  ;;
    450)  mu_range="2.0"  ;;
    500)  mu_range="2.0"  ;;
    550)  mu_range="1.0"  ;;
    600)  mu_range="0.8"  ;;
    700)  mu_range="0.6"  ;;
    800)  mu_range="0.6"  ;;
    900)  mu_range="0.6"  ;;
    1000) mu_range="0.6"  ;;
    1100) mu_range="0.6"  ;;
    1200) mu_range="0.6"  ;;
    1400) mu_range="0.6"  ;;
    1600) mu_range="0.6"  ;;
    *) { echo "Unknown mass: '${mass}'"; exit 1; } ;;
esac

# Set up PATH
script_dir="$(readlink -e bbtt_global_significance/scripts)"
PATH="${script_dir}:${PATH}"

outfile="${outdir}/toy_${seed}.csv"
runDiscoveryTestStatToys.py \
    "${infile}" \
    -s "${seed}" \
    -n "${ntoys}" \
    -o "${outfile}" \
    --mu-range "${mu_range}" \
    --optimizer-strategy 1 \
    2>&1

echo "Finished: $(date)"
