#!/bin/bash
str=$1

#Define paths and variables
inpath=/cndd3/dburrows/DATA/te/rna/PE.fastq/trimmed
outpath=/cndd3/dburrows/DATA/te/rna/PE.bam/
id_arr=($inpath/Sample_1848-GABA/ $inpath/Sample_4332-GABA/)


#Loop over each sample
for i in ${id_arr[@]}
do
    echo Running $i
    r1=$i/R1-merge_val_1.fq.gz
    r2=$i/R2-merge_val_2.fq.gz

    mkdir $outpath/$(basename $i)
    cd $outpath/$(basename $i)

    STAR --genomeDir /cndd2/jchien/iGenome/STAR_gencode_v37 --readFilesIn  $r1 $r2 --outSAMunmapped None --outFilterType BySJout --outSAMattributes All --outFilterMultimapNmax 100 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat --runThreadN 16 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --winAnchorMultimapNmax 200 --outMultimapperOrder Random --outSAMmultNmax -1    
    samtools index Aligned.sortedByCoord.out.bam
    cp $CODE3/te_ageing/workspace.sh $outpath/$(basename $i)/log.workspace
done

