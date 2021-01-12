#!/bin/bash

ENVNAME="input_utils_test_env"
conda install python=3.7 yaml pyyaml -ymp $ENVNAME

APATH="${ENVNAME}/etc/conda/activate.d"
DPATH="${ENVNAME}/etc/conda/deactivate.d"
mkdir -p ${APATH}
mkdir -p ${DPATH}

cp env_setup/a.welcome.sh env_setup/a.env_vars.sh $APATH
cp env_setup/d.env_vars.sh $APATH

