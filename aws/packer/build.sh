#!/bin/bash

# create a timestamp for the build
BUILD_TIMESTAMP=$(date +%Y%m%d-%H%M%S)

# export logging variables
export PACKER_LOG=1
export PACKER_LOG_PATH="logs/packer-detailed-$BUILD_TIMESTAMP.log"

# create logs directory if it doesn't exist
[ ! -d "logs" ] && mkdir -p "logs"

# build the AMI
packer build -var "build_timestamp=$BUILD_TIMESTAMP" -var-file=build-variables.json ami-rnaseq.json
