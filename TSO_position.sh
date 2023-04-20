#!/bin/bash

#Define paths
data_path="/cndd3/dburrows/DATA/te/rna/fastq.subset.TSO/"

#Define variables 
shopt -s extglob
inp_arr=($(ls $data_path/*merge*))
#TSO=AACGCAGAGTAC

#Loop through each dataset 
for i in ${inp_arr[@]}
do 
    echo 'TSO search running for ' $i
    targets=($(gunzip -c $i | grep $TSO))
    out_arr=()

    for t in ${targets[@]}
    do
        rest=${t#*$TSO}
        out_arr+=($(( ${#t} - ${#rest} - ${#TSO} )))
    done
    echo ${out_arr[@]} > $i-TSO-pos.txt
    echo 'TSO search complete for ' $i
done
echo 'Done'


gunzip -c 1848-GLU-R2-trim.fq.gz | cat | awk '{if(NR%4==2) print length($1)}' > read_length.txt