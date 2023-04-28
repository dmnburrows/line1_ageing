#!/bin/bash

str=$1

#Define paths and variables
datapath=$pe_rna_fa
outpath=/cndd3/dburrows/DATA/te/rna/PE.fastq/merged/
id=$(find $datapath/ -maxdepth 1 -name "*$str*")
id_arr=($id)

#Loop over each sample
for i in ${id_arr[@]}
do

    echo $i
    R1=$(ls $(find $i/fastq/ -maxdepth 1 -name "*R1*"))
    R2=$(ls $(find $i/fastq/ -maxdepth 1 -name "*R2*"))
    R1_arr=($R1)
    R2_arr=($R2)
    echo ${R1_arr[0]} ${R1_arr[1]} 
    echo ${R2_arr[0]} ${R2_arr[1]} 

    # # #concatenate lanes together for each read
    sample=$(basename $i) 
    mkdir $outpath/$sample
    cat  ${R1_arr[0]} ${R1_arr[1]}  > $outpath/$sample/R1-merge.fastq.qz
    cat  ${R2_arr[0]} ${R2_arr[1]}  > $outpath/$sample/R2-merge.fastq.qz
done

