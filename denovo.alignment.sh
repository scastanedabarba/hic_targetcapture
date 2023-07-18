#!/bin/bash

source /usr/modules/init/bash

#import necessary modules
module load samtools
module load bwa
module load samblaster

FILE=$1

base=$(basename $FILE _R1.fastq.gz)

#specify path to input files
R1_in=../03_trimmed_reads/${base}_R1.fastq.gz
R2_in=../03_trimmed_reads/${base}_R2.fastq.gz

genome=../02_ref_genomes/ecoli_pb10.fasta 

bam_out=../05_bwa_alignments/denovo/unsorted/${base}_aligned.bam
sorted=../05_bwa_alignments/denovo/sorted/${base}_sorted_aligned.bam

bwa mem -t 32 -5SPY $genome $R1_in $R2_in | samblaster |samtools view -@ 32 -S -h -b -F 2304 > $bam_out
samtools sort -o $sorted -n -@ 32 $bam_out #sorted by name rather than position

echo $base completed