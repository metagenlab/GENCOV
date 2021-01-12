#!/bin/bash



printf "%s\n" \
    "Loaded environment: genbank to fasta and gff converter" \
    "Location: $CONDA_PREFIX"\
    " "\
    "$(conda list |& grep -E 'mummer|^art|^perl\s|^python|biopython|json|csv')" \
    " "\
    "This environment is used for genbank to gff conversion " \
    " "\
    "simply run:" \
    "       gb2fasta_gff --help"
