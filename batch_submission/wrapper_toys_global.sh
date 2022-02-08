#!/usr/bin/env bash
set -eu
(( $# == 3 )) || { echo "Usage: wrapper_toys_global.sh indir outdir nToy"; exit 1; }

echo "Start: $(date)"

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

tar -xzf /cephfs/user/s6crdeut/WSMaker_code_compiled.tar.gz \
    || { echo "Cannot get WSMaker code"; exit 1; }

# Build workspaces
(
    set +eu
    cd WSMaker_HH_bbtautau
    source setup.sh
    set -eu

    # Replace input and output dir (they are symlinked to cephfs)
    [[ -e inputs ]] && rm -r inputs
    [[ -e output ]] && rm -r output
    mkdir inputs output

    # Link input histograms (faster via CephFS)
    mkdir -p inputs/combined_inputs
    ln -s "${indir}"/*.root inputs/combined_inputs/

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

    for mass in 251 260 280 300 325 350 375 400 450 500 550 \
                    600 700 800 900 1000 1100 1200 1400 1600; do
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

        runDiscoveryTestStat.py \
            "WSMaker_HH_bbtautau/output/${ws_name}/workspaces/combined/${mass}.root" \
            -o /jwd/outputs/"toy_m${mass}.csv" \
            -m "${mass}" \
            -i "${nToy}" \
            --mu-range "${mu_range}" \
            --optimizer-strategy 1 \
            --globs-tree "WSMaker_HH_bbtautau/inputs/combined_inputs/toy_globs_${mass}.root" \
            --globs-index "${nToy}" \
            2>&1 > /dev/null \
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

echo "Finished: $(date)"
