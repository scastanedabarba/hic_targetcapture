#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda activate ~/miniconda3/envs/mamba_base/envs/extract_hosts

FILE=$1

FILE2=$(basename $FILE)

echo working on $FILE2

python3 enrichment_comparison.data_analysis.py $FILE2