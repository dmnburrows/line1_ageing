#!/bin/bash

#Define paths
ind_path="/cndd3/dburrows/DATA/te/rna/aligned.subset.no_TSO/"
out_path="/cndd3/dburrows/DATA/te/rna/L1EM_counts.subset.no_TSO"

#Define parameters
id=$1
cell=$2
echo 'Running: ' $id-$cell

#Run L1EM
mkdir $out_path/$id-$cell
cd $out_path/$id-$cell
bash -e ~/L1EM/run_L1EM.sh $ind_path/$id-$cell ~/L1EM /cndd3/dburrows/DATA/te/gtf/hg38.fa 