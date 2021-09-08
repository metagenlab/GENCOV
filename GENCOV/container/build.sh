#!/usr/bin/bash

# Build all Docker images and tag them to match description in snakemake rules

# hoelzerm@rki.de

for DOCKERFILE in */Dockerfile; do
    TOOL=$(grep 'ENV TOOL' $DOCKERFILE | awk '{print $3}')
    VERSION=$(grep 'ENV VERSION' $DOCKERFILE | awk '{print $3}')
    HASH=$(md5sum $DOCKERFILE | awk '{print $1}' | cut -c 1-7)
    if [ $TOOL ]; then
        cd $(dirname $DOCKERFILE)
        docker build -t rkibioinf/${TOOL}:${VERSION}--${HASH} .
        #--no-cache helps to build from scratch and to not share image slices that might be not available from the organization account (e.g. private images)
        cd ..
    fi
done