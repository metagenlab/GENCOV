#!/bin/bash
 
SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
TMPDIR="$(realpath $(mktemp --tmpdir='./' -d -t '.ci-tiny_clipping-XXXXXXXXXX'))"
CONDA_BASE=$(conda info --base)

echo "tiny clipping test:"
echo ""
echo "dir: ${TMPDIR}"
echo ""

echo "preparing workspace ..."
find ${SCRIPTPATH}/* -exec cp -r '{}' ${TMPDIR} \;
wget -q -O "${TMPDIR}/kraken_db.tar.gz" https://zenodo.org/record/3854856/files/GRCh38.p13_GBcovid19-2020-05-22.tar.gz?download=1
tar -C "${TMPDIR}" -zxvf "${TMPDIR}/kraken_db.tar.gz"
rm -Rf "${TMPDIR}/kraken_db.tar.gz"

ESCAPED_REPLACE=$(printf '%s\n' "${TMPDIR}" | sed -e 's/[\/&]/\\&/g')
sed -i "s/{TESTDIR}/${ESCAPED_REPLACE}/" "${TMPDIR}/test.config" "${TMPDIR}/sample.config.yaml"

eval "$(conda shell.bash hook)"
conda env create -q --force -f "${SCRIPTPATH}/../../envs/general.yaml" --prefix "${TMPDIR}/_env"
source activate "${TMPDIR}/_env/"

cd "${TMPDIR}"

echo "creating simple amplicon data ..."
"${SCRIPTPATH}/../../scripts/sim_amplicons.py" --pmin 21 --pmax 21 --cmin 50 --cmax 100 --len 150 --out test NC_045512.2.fasta
snakemake -s "${SCRIPTPATH}/../../ncov_minipipe.snake" --configfile test.config --use-conda --cores 1 --forcerun -p --notemp

echo "checking primer clipped reads ..."
if [[ $(zcat "${TMPDIR}/test_clipped_R1.fq.gz" | md5sum) = $(zcat "${TMPDIR}/Analyses/trimmed/tiny_test/tiny_test.R1.noprimer.fastq.gz" | md5sum) ]]; then
   echo "R1 clipping correct"
else
  echo "ERROR: R1 clipping not as expected" >&2
  exit 1
fi

if [[ $(zcat "${TMPDIR}/test_clipped_R2.fq.gz" | md5sum) = $(zcat "${TMPDIR}/Analyses/trimmed/tiny_test/tiny_test.R2.noprimer.fastq.gz" | md5sum) ]]; then
   echo "R2 clipping correct"
else
   echo "ERROR: R2 clipping not as expected" >&2
   exit 1
fi

echo "checking consensus sequence 01 ..."
"${SCRIPTPATH}/../../scripts/compare_seqs.py" --cutoff 1000 -d "${TMPDIR}/test_reference.fasta" "${TMPDIR}/Analyses/consensus/ambiguous/tiny_test.consensus.fasta"

echo "checking consensus sequence 02 ..."
"${SCRIPTPATH}/../../scripts/compare_seqs.py" --cutoff 1000 -d "${TMPDIR}/test_reference.fasta" "${TMPDIR}/Analyses/consensus/masked/tiny_test.masked_consensus.fasta"


cd ..





