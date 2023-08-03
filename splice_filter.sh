#!/bin/bash
samtools view -h 5pfilt-tss.bam | awk '/^@/ {print;next} !($6 ~ /N/)' | samtools view -bo 5pfilt-tss_nosplice.bam
