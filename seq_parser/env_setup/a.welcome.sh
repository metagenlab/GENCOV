#!/bin/bash



printf "%s\n" \
    "Loaded environment: test_input_utils" \
    "Location: $CONDA_PREFIX"\
    " "\
    "$(conda list |& grep -E '^bracken|^kraken|^perl\s|^python|yaml|json|csv')" \
    " "\
    "This environment is used for testing the input_utils.py module" \
    "fileparser.py is an example of how to use it"\
