#!/bin/bash

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
#TMPDIR="$(realpath $(mktemp --tmpdir='./' -d -t '.ci-tiny_clipping-XXXXXXXXXX'))"
TMPDIR="/home/stephan/Documents/tmp/.ci-tiny_clipping-WDTMS5kKUM"
CONDA_BASE=$(conda info --base)

echo "tiny clipping test:"
echo ""
echo "dir: ${TMPDIR}"
echo ""


#cp -R "${SCRIPTPATH}/." "${TMPDIR}/"
#wget -O "${TMPDIR}/kraken_db.tar.gz" https://zenodo.org/record/3854856/files/GRCh38.p13_GBcovid19-2020-05-22.tar.gz?download=1
#tar -C "${TMPDIR}" -zxvf "${TMPDIR}/kraken_db.tar.gz"
#rm -Rf "${TMPDIR}/kraken_db.tar.gz"

#ESCAPED_REPLACE=$(printf '%s\n' "${TMPDIR}" | sed -e 's/[\/&]/\\&/g')
#sed -i "s/{TESTDIR}/"${ESCAPED_REPLACE}"/" "${TMPDIR}/test.config" "${TMPDIR}/sample.config.yaml"

#eval "$(conda shell.bash hook)"
#conda env create --force -f "${SCRIPTPATH}/../../envs/general.yaml" --prefix "${TMPDIR}/_env"
#conda activate "${TMPDIR}/_env/"

cd "${TMPDIR}"

"${SCRIPTPATH}"/../../scripts/sim_amplicons.py --out test NC_045512.2.fasta
snakemake -s "${SCRIPTPATH}/../../ncov_minipipe.snake" --configfile test.config --use-conda --cores 1 --forcerun -p

cd ..





