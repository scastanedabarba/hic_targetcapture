#!/bin/bash

#SBATCH --job-name=subsample_align
#SBATCH -o subsample_align.out
#SBATCH -e subsample_align.error
#SBATCH -C ceph
#SBATCH -p reg
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=cast9836@vandals.uidaho.edu
#SBATCH --array=0-124

FILES=(../07_subsampling/subsampled/*R1.fastq)

srun bash math_model.align_subsampled.sh ${FILES[SLURM_ARRAY_TASK_ID]}