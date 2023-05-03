#!/bin/bash

#Define paths and variables
inpath=$DATA3/te/rna/PE.fastq/merged/
outpath=/cndd3/dburrows/DATA/te/rna/PE.fastq/trimmed
id_arr=($inpath/Sample_5077-GLU/ $inpath/Sample_4337-GABA/ $inpath/Sample_1823-GLU/ $inpath/Sample_5293-GLU/ $inpath/Sample_5184-GLU/)

#Loop over each sample
for i in ${id_arr[@]}
do
    echo $i
    read_a=($(ls $i))
    $CODE3/te_ageing/trim_galore.sh 8 3 20 0.1 $i/${read_a[0]} $i/${read_a[1]} $outpath/$(basename $i)
    cp $CODE3/te_ageing/workspace-1.sh $outpath/$(basename $i)/log.workspace
    cp $CODE3/te_ageing/trim_galore.sh $outpath/$(basename $i)/log.trim_galore

done

