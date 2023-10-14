#!/bin/bash

#Define paths
in=($(find /cndd/dburrows/DATA/public_data/ATEM_benchmark/ -maxdepth 1 -name "*EN*"))
out=/cndd/dburrows/DATA/public_data/TET_benchmark/

# loop over all files in the input directory
for i in ${in[@]}
do
    echo 'Running: ' $i
    base=$(basename $i)
    TEcount --sortByPos --format BAM --mode multi --stranded reverse -b $i/Aligned.sortedByCoord.out.bam??? --GTF /cndd3/dburrows/DATA/te/annotations/gencode/gencode.v37.annotation.hg38.gtf --TE /cndd3/dburrows/DATA/te/annotations/rmsk/rmsk.hg38.gtf --outdir $out/TET-$base???
done

echo 'Done'
