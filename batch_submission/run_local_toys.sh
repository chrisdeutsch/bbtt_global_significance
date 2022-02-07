#!/usr/bin/env bash
set -eu
for mass in 251 260 280 300 325 350 375 400 \
                450 500 550 600 700 800 900 \
                1000 1100 1200 1400 1600; do
    outdir="/cephfs/user/s6crdeut/q0_toys/2022_02_07_combined/comb_${mass}"

    echo "Creating output directory: ${outdir}"
    mkdir -p "${outdir}"

    echo "Submitting m=${mass} ($(date))..."
    cluster=$(condor_submit Mass="${mass}" SeedOffset=0 submission_toys_local.jdl \
                  | tail -n 1 \
                  | sed "s/^.*cluster //; s/.$//")
    echo "Submitted to cluster '${cluster}'"

    echo "Waiting for job to finish..."
    condor_wait "logs/log.${cluster}" || { echo "Job failed..."; exit 1; }

    sleep 60

    # Merge csv files
    echo "Merging toys..."
    (
        cd "${outdir}"

        {
            head -n 1 toy_0.csv
            find -name "toy_*.csv" -exec awk "{if(NR>1)print}" {} \;
        } > "toys_combined_${mass}.csv"

        rm toy_*.csv
    )

    # Deleting logs
    rm logs/"out.${cluster}."*
    rm logs/"err.${cluster}."*
done
