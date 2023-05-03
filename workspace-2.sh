#!/bin/bash
str=$1

#Define paths and variables
inpath=$DATA3/te/rna/PE.fastq/merged/
outpath=/cndd3/dburrows/DATA/te/rna/PE.fastq/trimmed
id_arr=($(find $inpath/ -maxdepth 1 -name "*Sample*$str*"))

#Loop over each sample
for i in ${id_arr[@]:20:${#id_arr[@]}}}
do
    echo $i
    read_a=($(ls $i))
    $CODE3/te_ageing/trim_galore.sh 8 3 20 0.1 $i/${read_a[0]} $i/${read_a[1]} $outpath/$(basename $i)
    cp $CODE3/te_ageing/workspace-2.sh $outpath/$(basename $i)/log.workspace
    cp $CODE3/te_ageing/trim_galore.sh $outpath/$(basename $i)/log.trim_galore

done