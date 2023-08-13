#!/bin/bash

#get name of file currently being worked on
library=$1

echo $library
base=$(basename $library _R1.fastq.gz)
echo $base

#specify paths to input and output files
genome=../02_ref_genomes/combined_genomes.fasta #bwa indexed reference genomes
R1_in=${base}_R1.fastq.gz #R1
R2_in=${base}_R2.fastq.gz #R2
R1_trimmed=../03_trimmed_reads/${base}_R1.fastq.gz #trimmed R1
R2_trimmed=../03_trimmed_reads/${base}_R2.fastq.gz #trimmed R2
bam_out=../05_bwa_alignments/enrichment_comparison/unsorted/${base}_aligned.bam #alignment file
sorted=../05_bwa_alignments/enrichment_comparison/sorted/${base}_sorted_aligned.bam #sorted alignment file (by name)

#run fastp, bwa, and samtools
fastp -i $R1_in -I $R2_in -o $R1_trimmed -O $R2_trimmed --thread 16 -l 50
bwa mem -t 32 -5SPY $genome $R1_trimmed $R2_trimmed | samblaster |samtools view -@ 32 -S -h -b -F 2316 > $bam_out
samtools sort -o $sorted -n -@ 32 $bam_out #sorted by name rather than position