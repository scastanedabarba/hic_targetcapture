#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda activate ~/miniconda3/envs/mamba_base/envs/extract_hosts

FILE=$1

base=$(basename $FILE _aligned_sorted.bam)

echo working on $FILE with base $base

#specify file paths
input=../05_bwa_alignment/de_novo/${base}_aligned_sorted.bam
read_df=../06_output/de_novo/${base}_reads.csv
pothosts=../06_output/de_novo/${base}_pothosts.fasta
kraken_report=../06_output/de_novo/${base}.kreport2
kraken_out=../06_output/de_novo/${base}_out.kraken2
taxonomy_out=../06_output/de_novo/${base}_lineage.txt

#run python script to generate csv files
python3 extract_hosts.py $input $read_df $pothosts

#run kraken and then output report and file with taxids
kraken2 --db ~/kraken_standard --threads 32 --confidence .5 --report $kraken_report $pothosts > $kraken_out

#use file with taxid to generate file with 
awk -F '\t' '{print $3}' $kraken_out \
    | taxonkit reformat -I 1 \
    | csvtk -H -t sep -f 2 -s ';' -R \
    | csvtk add-header -t -n taxid,kingdom,phylum,class,order,family,genus,species > $taxonomy_out
