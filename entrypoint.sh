#!/bin/bash
# Read in the file of environment settings
source activate covpipe_environment/
# Then run the CMD
exec "$@"