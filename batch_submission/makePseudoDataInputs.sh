#!/usr/bin/env bash
input_histdir="/cephfs/user/s6crdeut/Workspaces/2021_08_26_CONF/hists"
pseudodata_histdir="/cephfs/user/s6crdeut/WSMakerPseudoData/pd_only"
outdir="/cephfs/user/s6crdeut/WSMakerPseudoData/merged"

for histpath in "${input_histdir}"/*.root; do
    histname=$(basename "${histpath}")

    nominal_input="${histpath}"
    pseudodata_input="${pseudodata_histdir}/${histname}"
    outfile="${outdir}/${histname}"

    if [[ ! -e "${pseudodata_input}" ]]; then
        echo "Histogram ${histname} does not exist"
        continue
    fi

    [[ -e "${outfile}" ]] && rm "${outfile}"
    hadd "${outfile}" "${nominal_input}" "${pseudodata_input}" || { echo "Error hadd'ing"; exit 1; }
done
