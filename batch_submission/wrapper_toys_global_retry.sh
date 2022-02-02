#!/usr/bin/env bash
set -eu
(( $# == 4 )) || { echo "Usage: wrapper_toys_global.sh indir outdir nToy mass"; exit 1; }

echo "Start: $(date)"

indir="$1"
outdir="$2"
nToy="$3"
mass="$4"

[[ -d "${indir}" ]] || { echo "Indir does not exist"; exit 1; }
[[ -d "${outdir}" ]] || { echo "Outdir does not exist"; exit 1; }

mkdir -p /jwd/run
mkdir -p /jwd/outputs

cd /jwd/run

# Get code
tar -xzf /cephfs/user/s6crdeut/bbtt_global_significance.tar.gz \
    || { echo "Cannot get bbtt_global_significance"; exit 1; }

tar -xzf /cephfs/user/s6crdeut/WSMaker_code.tar.gz \
    || { echo "Cannot get WSMaker code"; exit 1; }

# Build workspaces
(
    set +eu
    cd WSMaker_HH_bbtautau
    source setup.sh
    rm -r build/*
    (cd build/ && cmake .. && make) || { echo "Failure compiling WSMaker"; exit 1; }
    set -eu

    # Replace input and output dir (they are symlinked to cephfs)
    [[ -e inputs ]] && rm -r inputs
    [[ -e output ]] && rm -r output
    mkdir inputs outputs

    # Copy input histograms
    mkdir -p inputs/combined_inputs
    cp "${indir}"/*.root inputs/combined_inputs/

    buildWorkspace4HH.py "combined_inputs" "combined_pseudodata${nToy}_m${mass}" \
                         --signal 2HDM --mass "${mass}" \
                         --pseudo-data "${nToy}" \
                         --override-regions \
                         "13TeV_TauHH_2tag2pjet_0ptv_LL_OS_PNN${mass}" \
                         "13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_${mass}" \
                         "13TeV_TauLHLTT_2tag2pjet_0ptv_2HDM_PNN_${mass}" \
                         "13TeV_TwoLepton_2tag2pjet_0ptv_ZllbbCR_mLL" \
        || { echo "Error building workspace"; exit 1; }
)

# Get q0
(
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

    # Set up PATH
    script_dir="$(readlink -e bbtt_global_significance/scripts)"
    PATH="${script_dir}:${PATH}"

    ws_name="combined_inputs.combined_pseudodata${nToy}_m${mass}"

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


    postfix="idx${nToy}_m${mass}.csv"

    # Retry Minuit2 strategy 2
    runDiscoveryTestStat.py \
        "WSMaker_HH_bbtautau/output/${ws_name}/workspaces/combined/${mass}.root" \
        -o /jwd/outputs/"retry_minuit2strat2_${postfix}" \
        -m "${mass}" \
        -i "${nToy}" \
        --mu-range "${mu_range}" \
        --optimizer-strategy 2 \
        --globs-tree "WSMaker_HH_bbtautau/inputs/combined_inputs/toy_globs_${mass}.root" \
        --globs-index "${nToy}" \
        2>&1 > /dev/null \
        || { echo "Error fitting toys"; exit 1; }

    # Retry Minuit strategy 1
    runDiscoveryTestStat.py \
        "WSMaker_HH_bbtautau/output/${ws_name}/workspaces/combined/${mass}.root" \
        -o /jwd/outputs/"retry_minuitstrat1_${postfix}" \
        -m "${mass}" \
        -i "${nToy}" \
        --mu-range "${mu_range}" \
        --optimizer Minuit \
        --optimizer-strategy 1 \
        --globs-tree "WSMaker_HH_bbtautau/inputs/combined_inputs/toy_globs_${mass}.root" \
        --globs-index "${nToy}" \
        2>&1 > /dev/null \
        || { echo "Error fitting toys"; exit 1; }

    # Retry Minuit strategy 2
    runDiscoveryTestStat.py \
        "WSMaker_HH_bbtautau/output/${ws_name}/workspaces/combined/${mass}.root" \
        -o /jwd/outputs/"retry_minuitstrat2_${postfix}" \
        -m "${mass}" \
        -i "${nToy}" \
        --mu-range "${mu_range}" \
        --optimizer Minuit \
        --optimizer-strategy 2 \
        --globs-tree "WSMaker_HH_bbtautau/inputs/combined_inputs/toy_globs_${mass}.root" \
        --globs-index "${nToy}" \
        2>&1 > /dev/null \
        || { echo "Error fitting toys"; exit 1; }


    # Try with mu range altered
    mu_range_doubled=$(bc <<< "2.0 * ${mu_range}")
    mu_range_halfed=$(bc <<< "0.5 * ${mu_range}")

    # Retry Minuit2 strategy 1 doubled mu-range
    runDiscoveryTestStat.py \
        "WSMaker_HH_bbtautau/output/${ws_name}/workspaces/combined/${mass}.root" \
        -o /jwd/outputs/"retry_minuit2strat1mu2_${postfix}" \
        -m "${mass}" \
        -i "${nToy}" \
        --mu-range "${mu_range_doubled}" \
        --optimizer Minuit2 \
        --optimizer-strategy 1 \
        --globs-tree "WSMaker_HH_bbtautau/inputs/combined_inputs/toy_globs_${mass}.root" \
        --globs-index "${nToy}" \
        2>&1 > /dev/null \
        || { echo "Error fitting toys"; exit 1; }

    # Retry Minuit2 strategy 1 halfed mu-range
    runDiscoveryTestStat.py \
        "WSMaker_HH_bbtautau/output/${ws_name}/workspaces/combined/${mass}.root" \
        -o /jwd/outputs/"retry_minuit2strat1mu0p5_${postfix}" \
        -m "${mass}" \
        -i "${nToy}" \
        --mu-range "${mu_range_halfed}" \
        --optimizer Minuit2 \
        --optimizer-strategy 1 \
        --globs-tree "WSMaker_HH_bbtautau/inputs/combined_inputs/toy_globs_${mass}.root" \
        --globs-index "${nToy}" \
        2>&1 > /dev/null \
        || { echo "Error fitting toys"; exit 1; }

    # Retry Minuit strategy 1 doubled mu-range
    runDiscoveryTestStat.py \
        "WSMaker_HH_bbtautau/output/${ws_name}/workspaces/combined/${mass}.root" \
        -o /jwd/outputs/"retry_minuitstrat1mu2_${postfix}" \
        -m "${mass}" \
        -i "${nToy}" \
        --mu-range "${mu_range_doubled}" \
        --optimizer Minuit \
        --optimizer-strategy 1 \
        --globs-tree "WSMaker_HH_bbtautau/inputs/combined_inputs/toy_globs_${mass}.root" \
        --globs-index "${nToy}" \
        2>&1 > /dev/null \
        || { echo "Error fitting toys"; exit 1; }

    # Retry Minuit strategy 1 halfed mu-range
    runDiscoveryTestStat.py \
        "WSMaker_HH_bbtautau/output/${ws_name}/workspaces/combined/${mass}.root" \
        -o /jwd/outputs/"retry_minuitstrat1mu0p5_${postfix}" \
        -m "${mass}" \
        -i "${nToy}" \
        --mu-range "${mu_range_halfed}" \
        --optimizer Minuit \
        --optimizer-strategy 1 \
        --globs-tree "WSMaker_HH_bbtautau/inputs/combined_inputs/toy_globs_${mass}.root" \
        --globs-index "${nToy}" \
        2>&1 > /dev/null \
        || { echo "Error fitting toys"; exit 1; }
)

cd /jwd/outputs/
tar -czf retry_idx${nToy}_m${mass}.tar.gz *.csv
cp retry_idx${nToy}_m${mass}.tar.gz "${outdir}/"

echo "Finished: $(date)"
