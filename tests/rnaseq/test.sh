#!/bin/bash
set -e

# activate conda environment
source /usr/local/env/conda/bin/activate

# display the current conda environment
echo "Current Conda environment:"
echo $CONDA_DEFAULT_ENV

conda info --envs

# run test scripts
Rscript library.R
echo "Running test script..."
