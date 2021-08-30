#!/usr/bin/env bash
for mass in 251 260 280 300 325 350 375 400 \
                450 500 550 600 700 800 900 \
                1000 1100 1200 1400 1600; do

    echo "Submitting m=${mass}..."
    cluster=$(condor_submit Mass="${mass}" SeedOffset=0 submission_toys_local.jdl \
                  | tail -n 1 \
                  | sed "s/^.*cluster //; s/.$//")
    echo "Submitted to cluster '${cluster}'"

    echo "Waiting for job to finish..."
    condor_wait "logs/log.${cluster}" || { echo "Job failed..."; exit 1; }

    # Merge csv files
    echo "Merging toys..."
    (
        cd "/cephfs/user/s6crdeut/WSMakerToys_q0/combined_${mass}"

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
