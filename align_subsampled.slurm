#!/bin/bash

#SBATCH --job-name=subsample_align
#SBATCH -o out.out
#SBATCH -e error.out
#SBATCH -p reg
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=cast9836@vandals.uidaho.edu
#SBATCH --array=0-124

FILES=(../03.5_subsamples_reads/subsampled/*R1.fastq)

srun bash align_subsampled.sh ${FILES[SLURM_ARRAY_TASK_ID]}