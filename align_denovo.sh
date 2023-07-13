#!/bin/bash

source /usr/modules/init/bash

#import necessary modules
module load samtools
module load bwa

FILE=$1

base=$(basename $FILE _R1.fastq.gz)

#specify path to input files
R1_in=../03_trimmed_read/${base}_R1.fastq.gz
R2_in=../03_trimmed_read/${base}_R2.fastq.gz

genome=../04_ref_genomes/ecoli_pb10.fasta 

bam_out=../05_bwa_alignment/de_novo/${base}_aligned_sorted.bam

bwa mem -t 8 -5SPY $genome $R1_in $R2_in |samtools view -@ 8 -S -h -b -F 2304 | samtools sort -n -o $bam_out

echo $base completed