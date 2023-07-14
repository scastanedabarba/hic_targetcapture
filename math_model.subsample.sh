#!/bin/bash

#input is a fastq file, subsamples this file and its pair at different sequencing depths for 5 trials and outputs the file with the trial and depth in the name

source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

FILE=$1

base=$(basename $FILE _R1.fastq)
R1=../07_subsampling/${base}_R1.fastq
R2=../07_subsampling/${base}_R2.fastq
echo $R1
echo $R2

for depth in 50000000 60000000 70000000 80000000 90000000
do
    echo subsampling library $FILE at depth $depth
    for trial in {1..5}
    do
        echo trial $trial started
        seed=$RANDOM
        outR1=../07_subsampling/subsampled/${base}_${depth}_${trial}_R1.fastq
        outR2=../07_subsampling/subsampled/${base}_${depth}_${trial}_R2.fastq
        ~/miniconda3/envs/mamba_base/bin/seqtk sample -2 -s $seed $R1 $depth > $outR1
        ~/miniconda3/envs/mamba_base/bin/seqtk sample -2 -s $seed $R2 $depth > $outR2
        echo trial $trial complete
    echo subsampling at $depth complete
    done
done
