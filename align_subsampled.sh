#!/bin/bash

source /usr/modules/init/bash

#import necessary modules
module load samtools
module load bwa

FILE=$1

base=$(basename $FILE _R1.fastq)

#specify path to input files
R1_in=../03_trimmed_read/subsampled/${base}_R1.fastq
R2_in=../03_trimmed_read/subsampled/${base}_R2.fastq

genome=../04_ref_genomes/combined_genomes.fasta 

bam_out=../05_bwa_alignment/subsampled/${base}_aligned_sorted.bam

bwa mem -t 8 -5SPY $genome $R1_in $R2_in |samtools view -@ 8 -S -h -b -F 2316 | samtools sort -n -o $bam_out

echo $base completed