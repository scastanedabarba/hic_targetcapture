#!/bin/bash

#SBATCH --job-name=denovo_align
#SBATCH -o denovo_alignment_%A_%a.out
#SBATCH -e denovo_alignment_%A_%a.error
#SBATCH -p reg
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=cast9836@vandals.uidaho.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-12

FILES=(../03_trimmed_reads/*R1.fastq.gz)

srun bash denovo.alignment.sh ${FILES[SLURM_ARRAY_TASK_ID]}head