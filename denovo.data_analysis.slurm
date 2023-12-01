#!/bin/bash

#SBATCH --job-name=denovo_analysis
#SBATCH -o denovo_analysis_%A-%a.out
#SBATCH -e denovo_analysis_%A-%a.error
#SBATCH -C ceph
#SBATCH -p reg
#SBATCH --cpus-per-task=2
#SBATCH --mem=75G
#SBATCH --mail-user=cast9836@vandals.uidaho.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-11

FILES=(../05_bwa_alignments/denovo/sorted/*_sorted_aligned.bam)

srun bash denovo.data_analysis.sh ${FILES[SLURM_ARRAY_TASK_ID]}