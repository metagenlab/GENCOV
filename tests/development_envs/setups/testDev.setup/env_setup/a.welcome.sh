#!/bin/bash



printf "%s\n" \
    "Loaded environment: dev - ncov-minipipe testing" \
    "Location: $CONDA_PREFIX"\
    " "\
    "$(conda list |& grep -E '^bracken|^art|^perl\s|^python|yaml|json|csv|vcftools|bedtools')" \
    " "\
    "This environment is used for developing a test suite for the " \
    "ncov-minipipe"\
