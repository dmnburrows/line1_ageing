#!/bin/bash
#This script runs trim galore to filter out low quality reads, remove adapter sequences and clip certain ends -> change parameter as needed! 
#See https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md#version-064
#See https://cutadapt.readthedocs.io/en/stable/guide.html

#Define variables
cores=$1 #number of cores
clip=$2 #number of bases to clip from 5' end of read 2
qc=$3 #quality score cutoff
err=$4 #maximum error rate
r1=$5 #read 1 fastq file input
r2=$6 #read 2 fastq file input
out=$7 

trim_galore --cores $cores --paired --clip_R2 $clip -q $qc -e $err --fastqc_args "-noextract" $r1 $r2 -o $out 



