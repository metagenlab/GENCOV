#!/bin/bash

set -eEuo pipefail

ENVNAME="test_dev_env"
SOURCE="$( realpath "${BASH_SOURCE[0]}" )"
SOURCE_DIR="$( dirname "$SOURCE" )"
ENV_SETUP="${SOURCE_DIR}/env_setup"
OUTDIR=${1-.}
ENVPATH="${OUTDIR}/${ENVNAME}"
APATH="${ENVPATH}/etc/conda/activate.d"
DPATH="${ENVPATH}/etc/conda/deactivate.d"


stop_setup(){
    unset ENVNAME SOURCE SOURCE_DIR ENV_SETUP OUTDIR ENVPATH APATH DPATH
    exit 0
}

conda install \
    -c defaults -c conda-forge -c bioconda \
    python=3.7 \
    yaml \
    rsync \
    pandas \
    fastq-tools \
    pyyaml \
    tabix \
    biopython \
    art \
    bbmap \
    bcftools \
    bedtools \
    -ymp "${ENVPATH}"

mkdir -p "${APATH}"
mkdir -p "${DPATH}"

cp "${ENV_SETUP}/a.welcome.sh" "${ENV_SETUP}/a.env_vars.sh" "${APATH}"
cp "${ENV_SETUP}/d.env_vars.sh" "${DPATH}"

# Customization

stop_setup

