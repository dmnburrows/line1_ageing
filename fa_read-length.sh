#!/bin/bash

in=$1 #fastq file
out=$2 #output file

#This script calculates the length of all reads in a fastq file
gunzip -c $in | cat | awk '{if(NR%4==2) print length($1)}' > $out
