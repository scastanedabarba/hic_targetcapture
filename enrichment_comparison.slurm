#!/bin/bash

#SBATCH --job-name=enrichment_comparison
#SBATCH -o array_test_%A_%a.out
#SBATCH -e array_test_%A_%a.error
#SBATCH -p eight
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=cast9836@vandals.uidaho.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-5

source /usr/modules/init/bash
module load fastqc
module load python/3.6.7

#change to directory where raw data is stored 
cd ../01_raw_data

#make directories for output files
mkdir ../03_trimmed_reads
mkdir -p ../04_fastqc/trimmed ../03_fastqc/untrimmed
mkdir -p ../05_bwa_alignments/enrichment_comparison/unsorted ../05_bwa_alignments/enrichment_comparison/sorted
mkdir -p ../06_ouput/enrichment_comparison

FILES=(*R1.fastq.gz)

#run array job for each file
srun bash enrichment_comparison.alignment.sh ${FILES[SLURM_ARRAY_TASK_ID]}

#run fastqc and multiqc
fastqc *fastq.gz -o ../04_fastqc/untrimmed/
fastqc ../03_trimmed_reads/*fastq.gz -o ../04_fastqc/trimmed/

multiqc ../04_fastqc/trimmed/* -n trimmed
multiqc ../04_fastqc/untrimmed/* -n untrimmed

#run python script for doing count data
python3 enrichment_comparison.data_analysis.py
