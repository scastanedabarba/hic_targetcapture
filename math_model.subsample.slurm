#!/bin/bash

#SBATCH --job-name=subsample
#SBATCH -o subsample_%A-%a.out
#SBATCH -e subsample_%A-%a.error
#SBATCH -C ceph
#SBATCH -p reg
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=cast9836@vandals.uidaho.edu
#SBATCH --array=0-4

FILES=(../07_subsampling/*2.fastq)

srun bash math_model.subsample.sh ${FILES[SLURM_ARRAY_TASK_ID]}