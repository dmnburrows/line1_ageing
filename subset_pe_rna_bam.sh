#!/bin/bash

#Define paths
rna_path="/datasets/Public_Datasets/Dracheva_PsychEncode_development/processed/rna_seq/"
py_path="/cndd3/dburrows/DATA/te/rna/aligned.subset.no_TSO/"

#Define variables
cell_arr=('GLU', 'GABA')
id_arr=($(tail -n 10 $py_path/male_GLU_df.csv | cut -d ',' -f2))

#Loop through and send over each file
for i in ${id_arr[@]}
do 
	for e in ${cell_arr[@]}
	do 
		curr=$(find $rna_path -maxdepth 1 -name "*$i*$e*")
		ln -s "$curr/Aligned.sortedByCoord.out.bam" $py_path/$i-$e
		samtools index $py_path/$i-$e
		echo 'Sym linked and indexed ' $i-$e

	done
done

echo 'Done'