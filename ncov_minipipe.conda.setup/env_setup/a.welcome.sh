#!/bin/bash



printf "%s\n" \
    "Loaded environment: ncov_minipipe a.k.a CovPipe" \
    "Location: $CONDA_PREFIX"\
    " "\
    "$(conda list |& grep -E 'snakemake|^yaml|^perl\s|^python|strictyaml|pandas|csv')" \
    " "\
    "This environment is used to control environment tools for the pipeline" \
    " "\
    "simply run:" \
    "       ncov_minipipe --help"
