#!/bin/bash

source /usr/modules/init/bash

#import necessary modules
module load samtools
module load bwa
module load samblaster

FILE=$1

base=$(basename $FILE _R1.fastq)

#specify path to input files
R1_in=../07_subsampling/subsampled/${base}_R1.fastq
R2_in=../07_subsampling/subsampled/${base}_R2.fastq

genome=../02_ref_genomes/combined_genomes.fasta 

bam_out=../05_bwa_alignments/math_model/unsorted/${base}_aligned.bam
sorted=../05_bwa_alignments/math_model/sorted/${base}_sorted_aligned.bam

bwa mem -t 32 -5SPY $genome $R1_in $R2_in | samblaster |samtools view -@ 32 -S -h -b -F 2316 > $bam_out
samtools sort -o $sorted -n -@ 32 $bam_out #sorted by name rather than position

echo $base completed