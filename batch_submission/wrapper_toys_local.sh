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

# Set range of mu
# Roughly 5 x maximum of hadhad toys
filename=$(basename "${infile}")
mass="${filename%.root}"

# Set sensible fit range determined from maximum fitted |muhat| in 500
# trials with factor 10 for safety (rounded up to the next tenth)
mu_range=""
case "${mass}" in
    251)  mu_range="8.0"  ;;
    260)  mu_range="15.8" ;;
    280)  mu_range="17.1" ;;
    300)  mu_range="12.6" ;;
    325)  mu_range="11.0" ;;
    350)  mu_range="7.1"  ;;
    375)  mu_range="4.0"  ;;
    400)  mu_range="3.2"  ;;
    450)  mu_range="1.2"  ;;
    500)  mu_range="1.1"  ;;
    550)  mu_range="0.7"  ;;
    600)  mu_range="0.6"  ;;
    700)  mu_range="0.4"  ;;
    800)  mu_range="0.3"  ;;
    900)  mu_range="0.3"  ;;
    1000) mu_range="0.2"  ;;
    1100) mu_range="0.4"  ;;
    1200) mu_range="0.3"  ;;
    1400) mu_range="0.4"  ;;
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
    --optimizer-strategy 0 \
    2>&1
