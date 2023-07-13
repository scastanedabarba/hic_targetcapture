#!/bin/bash

#import necessary modules
module load samtools
module load bwa
module load samblaster
module load fastqc
module load fastp

#move to directory where raw data is
cd ../01_raw_data

#
library=$1
base=$(basename $libary _R1.fastq.gz)

#specify paths to input and output files
genome=../02_ref_genomes/combined_genomes.fasta #bwa indexed reference genomes
R1_in=../01_raw_data/${base}_R1.fastq.gz #R1
R2_in=../01_raw_data/${base}_R2.fastq.gz #R2
R1_trimmed=../03_trimmed_reads/${base}_R1.fastq.gz #trimmed R1
R2_trimmed=../03_trimmed_reads/${base}_R2.fastq.gz #trimmed R2
bam_out=../05_bwa_alignments/enrichment_comparison/${base}_aligned.bam #alignment file
sorted=../05_bwa_alignments/enrichment_comparison/sorted/${base}_sorted_aligned.bam #sorted alignment file (by name)

#run fastp, bwa, and samtools
fastp -i $R1_in -I $R2_in -o $R1_trimmed -O $R2_trimmed --thread 16 -l 50
bwa mem -t 32 -5SPY $genome $R1_trimmed $R2_trimmed | samblaster |samtools view -@ 32 -S -h -b -F 2304 > $bam_out
samtools sort -o $sorted -n -@ 32 $bam_out #sorted by name rather than position