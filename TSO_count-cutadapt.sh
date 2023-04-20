#!/bin/bash


#Define paths
fa_path="/cndd2/Public_Datasets/Dracheva_PsychEncode_development/raw_May2022/"
out_path="/cndd3/dburrows/DATA/te/rna/fastq.subset.TSO/"
py_path="/cndd3/dburrows/DATA/te/rna/aligned.subset.no_TSO/"

#Define variables
id=$1
cell=$2
TSO=AACGCAGAGTAC

#TSO filter
echo 'TSO filtering ' $id-$cell
cutadapt -g $TSO -G $TSO --overlap 10 --pair-filter=both --discard-untrimmed -e 0.2 -o $out_path/$id-$cell-R1-trim-rel.fq.gz -p $out_path/$id-$cell-R2-trim-rel.fq.gz $out_path/$id-$cell-R1-merge.fq $out_path/$id-$cell-R2-merge.fq 
echo 'Done ' $id-$cell

