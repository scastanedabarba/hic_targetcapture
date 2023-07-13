#!/bin/bash

#SBATCH --job-name=subsample
#SBATCH -o out_%A-%a.out
#SBATCH -e error_%A-%a.out
#SBATCH -p reg
#SBATCH --mail-user=cast9836@vandals.uidaho.edu
#SBATCH --array=0-4

FILES=(../03.5_subsamples_reads/*R1.fastq)

srun bash subsampling.sh ${FILES[SLURM_ARRAY_TASK_ID]}