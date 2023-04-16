#!/bin/bash

#Define paths
fa_path="/cndd2/Public_Datasets/Dracheva_PsychEncode_development/raw_May2022/"
out_path="/cndd3/dburrows/DATA/te/rna/fastq.subset.TSO/"
py_path="/cndd3/dburrows/DATA/te/rna/aligned.subset.no_TSO/"

#Define variables
id_arr=($(tail -n 10 $py_path/male_GLU_df.csv | cut -d ',' -f2))
cell=$1

#Loop through and send over each file
for id in ${id_arr[@]}
do 
    #Define paths
    echo 'Running: ' $id-$cell

    #Symlink
    curr=$(find $fa_path -maxdepth 1 -name "*$id*$cell*")
    R1=$(ls $(find $curr/fastq/ -maxdepth 1 -name "*R1*"))
    R2=$(ls $(find $curr/fastq/ -maxdepth 1 -name "*R2*"))
    R1_arr=($R1)
    R2_arr=($R2)
    ln -s "${R1_arr[0]}" $out_path/$id-$cell-R1-L1.fq
    ln -s "${R1_arr[1]}" $out_path/$id-$cell-R1-L2.fq
    ln -s "${R2_arr[0]}" $out_path/$id-$cell-R2-L1.fq
    ln -s "${R2_arr[1]}" $out_path/$id-$cell-R2-L2.fq
    echo 'Sym linked ' $id-$cell

    #Cat
    cd 
    echo 'Concatenating ' $id-$cell
    cat  $out_path/$id-$cell-R1-L1.fq $out_path/$id-$cell-R1-L2.fq > $out_path/$id-$cell-R1-merge.fq
    cat  $out_path/$id-$cell-R2-L1.fq $out_path/$id-$cell-R2-L2.fq > $out_path/$id-$cell-R2-merge.fq

    #TSO filter
    echo 'TSO filtering ' $id-$cell
    cutadapt -g GGTATCAACGCAGAGTACG -G GGTATCAACGCAGAGTACG --overlap 10 --pair-filter=both --discard-untrimmed  -o $out_path/$id-$cell-R1-trim.fq.gz -p $out_path/$id-$cell-R2-trim.fq.gz $out_path/$id-$cell-R1-merge.fq $out_path/$id-$cell-R2-merge.fq > $out_path/$id-$cell-log.txt

    #STAR map
    echo 'aligning ' $id-$cell
    mkdir /cndd3/dburrows/DATA/te/rna/aligned.subset.TSO/$id-$cell
    cd /cndd3/dburrows/DATA/te/rna/aligned.subset.TSO/$id-$cell
    STAR --genomeDir /cndd3/dburrows/DATA/te/gtf/STAR.gencode.v37 --readFilesIn $out_path/$id-$cell-R1-trim.fq.gz $out_path/$id-$cell-R2-trim.fq.gz --outSAMunmapped None --outFilterType BySJout --outSAMattributes All --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat --runThreadN 16 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts

    #TE transcripts
    #echo 'TET' $id-$cell
    #conda activate tet
    #TEcount --sortByPos --format BAM --mode multi -b /cndd3/dburrows/DATA/te/rna/aligned.subset.TSO/$id-$cell/Aligned.sortedByCoord.out.bam --GTF /cndd3/dburrows/DATA/te/gtf/gencode.v37.annotation.hg38.gtf --TE /cndd3/dburrows/DATA/te/gtf/rmsk.hg38.gtf --project $id-$cell --outdir /cndd3/dburrows/DATA/te/rna/tet_counts.subset.TSO
done

echo 'Done'

