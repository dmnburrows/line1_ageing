#!/bin/bash

#Define paths
rna_path="/cndd2/Public_Datasets/Dracheva_PsychEncode_development/raw_May2022/"
py_path="/cndd3/dburrows/DATA/te/rna/fastq.subset.no_TSO/"

#Define variables
cell_arr=('GLU' 'GABA')
id_arr=($(tail -n 10 $py_path/male_GLU_df.csv | cut -d ',' -f2))

#Loop through and send over each file
for i in ${id_arr[@]}
do 
	for e in ${cell_arr[@]}
	do 
		curr=$(find $rna_path -maxdepth 1 -name "*$i*$e*")
		file=$(find $curr/fastq/ -maxdepth 1 -name "*R2*")
		file_arr=($file)
		first_file=${file_arr[0]}
		ln -s "$first_file" $py_path/$i-$e
		echo 'Sym linkedd ' $i-$e

	done
done

echo 'Done'