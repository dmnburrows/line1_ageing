#!/bin/bash
cell=$1
parallel /cndd3/dburrows/CODE/te_ageing/batch_L1EM.sh {1} $cell :::: /cndd3/dburrows/DATA/te/rna/aligned.subset.no_TSO/id_vec.txt 