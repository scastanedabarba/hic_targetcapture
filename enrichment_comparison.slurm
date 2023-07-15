#!/bin/bash

#SBATCH --job-name=enrichment_comparison
#SBATCH -o array_test_%A_%a.out
#SBATCH -e array_test_%A_%a.error
#SBATCH -p eight
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=cast9836@vandals.uidaho.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-12

source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

cd ../01_raw_data/
FILES=(*R1.fastq.gz)

#run array job for each file
srun bash ../00_scripts/enrichment_comparison.alignment.sh ${FILES[SLURM_ARRAY_TASK_ID]}

#run fastqc and multiqc
#fastqc *fastq.gz -o ../04_fastqc/untrimmed/
#fastqc ../03_trimmed_reads/*fastq.gz -o ../04_fastqc/trimmed/

#multiqc ../04_fastqc/trimmed/* -n trimmed
#multiqc ../04_fastqc/untrimmed/* -n untrimmed