#!/bin/bash

#SBATCH --job-name=subsample
#SBATCH -o out_%A-%a.out
#SBATCH -e error_%A-%a.out
#SBATCH -p reg
#SBATCH --mail-user=cast9836@vandals.uidaho.edu
#SBATCH --array=0-4

#make necessary directories
mkdir -p ../07_subsampling/subsampled
FILES=(../07_subsampling/*R1.fastq)

srun bash math_model.subsample.sh ${FILES[SLURM_ARRAY_TASK_ID]}