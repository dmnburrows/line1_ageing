#!/bin/bash

#Define paths
data_path="/cndd3/dburrows/DATA/te/rna/fastq.subset.TSO/"

#Define variables 
output=$1
cell_arr=('GLU' 'GABA')
id_arr=($(tail -n 10 /cndd3/dburrows/DATA/te/rna/aligned.subset.no_TSO/male_GLU_df.csv | cut -d ',' -f2))
echo "ID,counts, total reads, normalised counts" > $output

#Loop through each dataset and calculate TSO counts
#NAME, N COUNTS, NORMALISED COUNTS
for i in ${id_arr[@]}
do 
	for e in ${cell_arr[@]}
	do 
        echo 'TSO counts running for ' $i-$e
        counts=$(gunzip -c $data_path/$i-$e-R1-trim-rel.fq.gz | wc -l)
        total=$(gunzip -c $data_path/$i-$e-R1-merge.fq | wc -l)
        echo $i-$e, $counts, $total, $(($counts/$total)) >> $output
		echo 'TSO counts complete for ' $i-$e

	done
done

echo 'Done'

