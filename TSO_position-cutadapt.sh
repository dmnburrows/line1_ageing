#!/bin/bash

#Define paths
data_path="/cndd3/dburrows/DATA/te/rna/fastq.subset.TSO/"

#Define variables 
inp_arr=($(ls $data_path/*trim*rel*gz))

#Loop through each dataset 
for i in ${inp_arr[@]}
do 
    echo 'TSO search running for ' $i
    gunzip -c $i | cat | awk '{if(NR%4==2) print length($1)}' > $i-read_length.txt
    echo 'TSO search complete for ' $i
done
echo 'Done'
