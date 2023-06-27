#!/bin/bash

#Define paths
# clust=($(ls $in_path))

# loop over all files in the input directory
i=$1
in_path="/cndd3/dburrows/DATA/public_datasets/10x.NSCLC_tumour.5p/sep_bam/"
echo 'Running: ' $i
TEcount --format BAM --mode multi -b $in_path/$i/Aligned.sortedByCoord.out.bam --GTF /cndd3/dburrows/DATA/te/gtf/annotations/gencode/gencode.v37.annotation.hg38.gtf --TE /cndd3/dburrows/DATA/te/gtf/annotations/rmsk/rmsk.hg38.gtf --project TET-$i --outdir /cndd3/dburrows/DATA/public_datasets/10x.NSCLC_tumour.5p/TET/

echo 'Done'

