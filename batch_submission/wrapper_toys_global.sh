#!/usr/bin/env bash
(( $# == 3 )) || { echo "Usage: wrapper_toys_global.sh indir outdir nToy"; exit 1; }

indir="$1"
outdir="$2"
nToy="$3"

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
    cd WSMaker_HH_bbtautau
    source setup.sh
    rm -r build/*
    (cd build/ && cmake .. && make) || { echo "Failure compiling WSMaker"; exit 1; }

    # Replace input and output dir (they are symlinked to cephfs)
    rm -r inputs output
    mkdir inputs output

    # Copy input histograms
    mkdir -p inputs/combined_inputs
    cp "${indir}"/*.root inputs/combined_inputs/

    for mass in 251 260 280 300 325 350 375 400 450 500 550 \
                600 700 800 900 1000 1100 1200 1400 1600; do
        buildWorkspace4HH.py "combined_inputs" "combined_pseudodata${nToy}_m${mass}" \
                             --signal 2HDM --mass "${mass}" \
                             --pseudo-data "${nToy}" \
                             --override-regions \
                             "13TeV_TauHH_2tag2pjet_0ptv_LL_OS_PNN${mass}" \
                             "13TeV_TauLH_2tag2pjet_0ptv_2HDM_PNN_${mass}" \
                             "13TeV_TauLHLTT_2tag2pjet_0ptv_2HDM_PNN_${mass}" \
                             "13TeV_TwoLepton_2tag2pjet_0ptv_ZllbbCR_mLL" \
            || { echo "Error building workspace"; exit 1; }
    done
)

# Get q0
(
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

    for mass in 251 260 280 300 325 350 375 400 450 500 550 \
                    600 700 800 900 1000 1100 1200 1400 1600; do
        ws_name="combined_inputs.combined_pseudodata${nToy}_m${mass}"

        runDiscoveryTestStat.py \
            "WSMaker_HH_bbtautau/output/${ws_name}/workspaces/combined/${mass}.root" \
            -o /jwd/outputs/"toy_m${mass}.csv" \
            -m "${mass}" \
            -i "${nToy}" \
            || { echo "Error fitting toys"; exit 1; }
    done
)

# Merge CSV
cd /jwd/outputs
{
    head -n 1 toy_m251.csv
    find -name "toy_*.csv" -exec awk "{if(NR>1)print}" {} \;
} > "toys_${nToy}.csv"

cp "toys_${nToy}.csv" "${outdir}/"
