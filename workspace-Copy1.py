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
    echo $start
    samtools view -b -F 0x40 {in_path}/{d}/Aligned.sortedByCoord.out.bam > {out_path}/{d}/Aligned.sortedByCoord.out.r2.bam
    echo $samtools filter r2
    samtools index {out_path}/{d}/Aligned.sortedByCoord.out.r2.bam
    samtools view -b -F 0x100 {out_path}/{d}/Aligned.sortedByCoord.out.r2.bam > {out_path}/{d}/Aligned.sortedByCoord.out.r2.unq.bam
    echo $
    samtools index {out_path}/{d}/Aligned.sortedByCoord.out.r2.unq.bam
    bedtools coverage -s -a /cndd3/dburrows/DATA/annotations/rmsk/l1hs_tss_bind.bed -b {out_path}/{d}/Aligned.sortedByCoord.out.r2.unq.bam > {out_path}/{d}/l1hs_unq.cov 
    bedtools coverage -s -a /cndd3/dburrows/DATA/annotations/rmsk/l1hs_tss_bind.bed -b {out_path}/{d}/Aligned.sortedByCoord.out.r2.bam > {out_path}/{d}/l1hs.cov 
    """

    get_ipython().run_cell_magic('bash', '', run)
    print('Done')
    

#this section shows what default arguments will be run if just executing the script
if __name__ == "__main__":
    in_path = '/cndd3/dburrows/DATA/te/rna/PE.bam/'
    out_path ='/cndd/dburrows/DATA/te/rna/PE.counts/l1hs_count_bind/'
    dir_list = [os.path.basename(i) for i in glob.glob(in_path + '*Samp*')]
    
    # Create a pool of 10 worker processes
    with multiprocessing.Pool(1) as pool:
        # Map the process_directory function to each item in dir_list
        pool.map(process_directory, dir_list)

    print("All tasks completed.")
    