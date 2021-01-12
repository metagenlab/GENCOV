#!/bin/bash

NEW_PREFIX="tool-env: genbank2fasta"
export PS1="$(sed "s;""$CONDA_DEFAULT_ENV"";$NEW_PREFIX;g"  <<< $PS1)"
export CONDA_PROMPT_MODIFIER="($NEW_PREFIX) "
