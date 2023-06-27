#!/bin/bash

in_path="/cndd3/dburrows/DATA/public_datasets/10x.NSCLC_tumour.5p/sep_bam/"
curr_list=($(ls $in_path))
for i in ${curr_list[@]}
do  
    echo 'Running: ' $i
    samtools view -h  $in_path/$i/Aligned.sortedByCoord.out-old.bam | awk 'BEGIN{FS=OFS="\t"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0}  $3~/^[1-9]|X|Y|MT/{$3="chr"$3; print $0} ' | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | samtools view -bS - > $in_path/$i/Aligned.sortedByCoord.out.bam


done
echo 'Done'
