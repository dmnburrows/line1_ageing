#CONDA ENV base_conda (python 3.9.7)
#Import packages
#---------------------------------------
import sys
import os
import glob
import pandas as pd
import numpy as np

#use pool multiprocessing
import multiprocessing
import glob


def process_directory(d):
    run = f"""
    samtools view -b -F 0x100 {out_path}/{d}/L1HS_all.bam > {out_path}/{d}/L1HS_all.unq.bam
    samtools index {out_path}/{d}/L1HS_all.unq.bam
    bedtools intersect -s -u -a {in_path}/{d}/Aligned.sortedByCoord.out.r2.bam -b /cndd3/dburrows/DATA/annotations/rmsk/rmsk.hg38.filt-5ptrim.merge.L1HS_all.canon.extended.bed > {out_path}/{d}/L1HS_all.extended.bam
    samtools index {out_path}/{d}/L1HS_all.extended.bam
    samtools view -b -F 0x100 {out_path}/{d}/L1HS_all.extended.bam > {out_path}/{d}/L1HS_all.extended.unq.bam
    samtools index {out_path}/{d}/L1HS_all.extended.unq.bam

    """
    get_ipython().run_cell_magic('bash', '', run)


#this section shows what default arguments will be run if just executing the script
if __name__ == "__main__":
    in_path = '/cndd3/dburrows/DATA/te/rna/PE.bam/'
    out_path='/cndd/dburrows/DATA/te/rna/PE.counts/l1hs/'
    dir_list = [os.path.basename(i) for i in glob.glob(in_path + '*Samp*')]
    
    # Create a pool of 10 worker processes
    with multiprocessing.Pool(10) as pool:
        # Map the process_directory function to each item in dir_list
        pool.map(process_directory, dir_list)

    print("All tasks completed.")
    