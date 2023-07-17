#!/bin/bash

#SBATCH --job-name=identify_hosts
#SBATCH -o array_test_%A_%a.out
#SBATCH -e array_test_%A_%a.error
#SBATCH -p eight
#SBATCH --cpus-per-task=16
#SBATCH --mem=70G
#SBATCH --mail-user=cast9836@vandals.uidaho.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-5

FILES=(../05_bwa_alignment/de_novo/*+_aligned*)

srun bash identify_hosts.sh ${FILES[SLURM_ARRAY_TASK_ID]}