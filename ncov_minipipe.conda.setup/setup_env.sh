#!/bin/bash

set -eEuo pipefail

SKIPCONDAINSTALL=0
while getopts s flag; do
    case "$flag" in 
        s) SKIPCONDAINSTALL=1
           ;;
    esac  
done
shift $(( $OPTIND -1 )) 

ENVNAME="covpipe_environment"
SOURCE="$( realpath "${BASH_SOURCE[0]}" )"
SOURCE_DIR="$( dirname "$SOURCE" )"
ENV_SETUP="${SOURCE_DIR}/env_setup"
OUTDIR=${1-.}
ENVPATH="${OUTDIR}/${ENVNAME}"
APATH="${ENVPATH}/etc/conda/activate.d"
DPATH="${ENVPATH}/etc/conda/deactivate.d"
OPT_PATH="${ENVPATH}/opt/ncov_minipipe"
SNAKE_DIR="$( realpath "${SOURCE_DIR}/.." )"

echo $OUTDIR

stop_setup(){
    unset ENVNAME SOURCE SOURCE_DIR ENV_SETUP OUTDIR ENVPATH APATH DPATH
    unset OPT_PATH SNAKE_DIR SKIPCONDAINSTALL
    exit 0
}

cp_if_file_exists(){
	if [ -f "$1" ] || [ -d "$1" ]; then
		cp -fr "$1" ${2-""}  
	fi
}

update_conda_or_install(){
    conda install \
        -c conda-forge  -c bioconda -c defaults \
        "python>=3.6.0" \
        "snakemake>=5.26" \
        strictyaml \
        -ymp "${ENVPATH}"
}

[ $SKIPCONDAINSTALL -eq 1 ] || update_conda_or_install

mkdir -p "${APATH}"
mkdir -p "${DPATH}"

cp "${ENV_SETUP}/a.welcome.sh" "${ENV_SETUP}/a.env_vars.sh" "${APATH}"
cp "${ENV_SETUP}/d.env_vars.sh" "${DPATH}"

# Customization

mkdir -p "$OPT_PATH"
#cp -r "${SNAKE_DIR}/"{adapter.fasta,envs,ncov_minipipe.py,ncov_minipipe.Rmd,seq_parser,tests} "${OPT_PATH}/." #tests contains tests right now
cp -fr "${SNAKE_DIR}/"{envs,scripts,ncov_minipipe.py,__version__.py,ncov_minipipe.Rmd,ncov_minipipe.snake,ncov_minipipe.config,seq_parser,rules} "${OPT_PATH}/."
cp_if_file_exists .git "${OPT_PATH}/."
cp_if_file_exists adapter.fasta "${OPT_PATH}/."
ln -frs "${OPT_PATH}/ncov_minipipe.py" "${ENVPATH}/bin/ncov_minipipe"


stop_setup

