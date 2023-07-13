#!/bin/bash

#SBATCH --job-name=denovo_align
#SBATCH -o out.out
#SBATCH -e error.out
#SBATCH -p reg
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=cast9836@vandals.uidaho.edu
#SBATCH --array=0-12

FILES=(../03_trimmed_read/*R1.fastq.gz)

srun bash align_denovo.sh ${FILES[SLURM_ARRAY_TASK_ID]}