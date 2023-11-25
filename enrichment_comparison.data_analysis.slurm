#!/bin/bash

#SBATCH --job-name=data_analysis
#SBATCH -o data_analysis_%A-%a.out
#SBATCH -e data_analysis_%A-%a.error
#SBATCH -p reg
#SBATCH -C ceph
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=cast9836@vandals.uidaho.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-11

FILES=(../05_bwa_alignments/enrichment_comparison/sorted/*.bam)

srun bash enrichment_comparison.data_analysis.sh ${FILES[SLURM_ARRAY_TASK_ID]}