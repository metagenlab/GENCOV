#!/bin/bash

NEW_PREFIX="tool-env: ncov_minipipe"
export PS1="$(sed "s;""$CONDA_DEFAULT_ENV"";$NEW_PREFIX;g"  <<< $PS1)"
export CONDA_PROMPT_MODIFIER="($NEW_PREFIX) "
