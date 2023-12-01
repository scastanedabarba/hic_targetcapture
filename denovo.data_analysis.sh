#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda activate ~/miniconda3/envs/mamba_base/envs/extract_hosts

FILE=$1

FILE2=$(basename $FILE)

echo working on $FILE2

#run python script to generate csv files
python3 denovo.data_analysis.py $FILE2

#specify file paths for kraken ouput
PREFIX=$(basename $FILE2 _sorted_aligned.bam)
pothosts=../06_output/denovo/${PREFIX}_pothosts.fasta
kraken_report=../06_output/denovo/${PREFIX}.kreport2
kraken_out=../06_output/denovo/${PREFIX}_out.kraken2
taxonomy_out=../06_output/denovo/${PREFIX}_lineage.txt

#run kraken and then use taxonkit+csvtk to generate file with taxonomic breakdown
kraken2 --db ~/kraken_standard --threads 32 --confidence .5 --report $kraken_report $pothosts > $kraken_out

awk -F '\t' '{print $3}' $kraken_out \
    | taxonkit reformat -I 1 \
    | csvtk -H -t sep -f 2 -s ';' -R \
    | csvtk add-header -t -n taxid,kingdom,phylum,class,order,family,genus,species > $taxonomy_out