#!/bin/bash

#SBATCH --job-name=clean_and_align
#SBATCH -o stdout.out
#SBATCH -e stderr.out
#SBATCH -p reg
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cast9836@vandals.uidaho.edu

source /usr/modules/init/bash

#import necessary modules
module load samtools
module load bwa
module load samblaster
module load fastqc
module load python/3.6.7
module load fastp

#script assumes that there is already a 00_scripts, 01_raw_data, and 04_ref_genomes directory
#   00_scripts contains this script, this is the directory from which it runs
#   01_raw_data contains that raw paired end sequencing files
#   04_ref_genomes contains the already indexed references to be used for bwa

#steps included here are
#   1) trim and clean reads using fastp, these are stored in the 03_trimmed_reads directory
#   2) run bwa alignment on each pair of R1 and R2 files, run the aligned file through samblaster, then filter using samtools
#   3) sort reads using samtools
#   4) generate fastqc files on trimmed and untrimmed reads, generate multiqc reports as well


#run commands from within 01_raw_data directory
cd ../01_raw_data

#specify path to indexed reference genomes
genome=../04_ref_genomes/combined_genomes.fasta 

#loop over files in 01_raw_data
for file in ~/00_projects/hic_targetcapture/01_raw_data/*R1.fastq.gz
    do
    base=$(basename $file _R1.fastq.gz)
    echo "working with sample $base"

    #specify path to input files
    R1_in=../01_raw_data/${base}_R1.fastq.gz
    R2_in=../01_raw_data/${base}_R2.fastq.gz

    #trim reads
    #mkdir ../03_trimmed_read
    R1_trimmed=../03_trimmed_read/${base}_R1.fastq.gz
    R2_trimmed=../03_trimmed_read/${base}_R2.fastq.gz
    #fastp -i $R1_in -I $R2_in -o $R1_trimmed -O $R2_trimmed --thread 16 -l 50

    #bwa alignment
    #mkdir -p ../05_bwa_alignment
    bam_out=../05_bwa_alignment/${base}_aligned.bam
    bwa mem -t 32 -5SPY $genome $R1_trimmed $R2_trimmed | samblaster |samtools view -@ 32 -S -h -b -F 2304 > $bam_out

    #sort reads
    #mkdir ../05_bwa_alignment/sorted
    sorted=../05_bwa_alignment/sorted/${base}_sorted_aligned.bam
    samtools sort -o $sorted -n -@ 32 $bam_out #sorted by name rather than position
    done

#run fastqc on untrimmed reads
#mkdir -p ../02_fastqc/untrimmed
#fastqc ../01_raw_data/*fastq.gz -o ../02_fastqc/untrimmed/
#multiqc ../02_fastqc/untrimmed/* -n untrimmed

#run fastqc on trimmed reads
#mkdir -p ../02_fastqc/trimmed
#fastqc ../03_trimmed_read/*fastq.gz -o ../02_fastqc/trimmed/
#multiqc ../02_fastqc/trimmed/* -n trimmed