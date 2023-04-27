#!/bin/bash

#!/bin/bash
celltype=$1
parallel /cndd3/dburrows/CODE/te_ageing/TSO_count-cutadapt.sh {1} $celltype :::: /cndd3/dburrows/DATA/te/rna/aligned.subset.no_TSO/id_vec.txt 