#!/bin/bash
#This script runs cutadapt to only keep sequences that contain a sequence of interest

#Define variables
in_r1=$1 #read 1 fastq file input
in_r2=$2 #read 2 fastq file input
out_r1=$3 #read 1 filtered output
out_r2=$4 #read 2 filtered output
TSO=$5 #sequence of interest
ov=$6 # minimum overlap
err=$7 # maximum error rate


cutadapt -g $TSO -G $TSO --overlap $ov --pair-filter=both --discard-untrimmed -e $err -o $out_r1 -p $out_r2 $in_r1 $in_r2 
