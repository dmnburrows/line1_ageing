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
    # The commands to be executed for each directory
    run = f""" bedtools intersect -s -sorted -wb -b {bam_path}/{d}/Aligned.sortedByCoord.out.bam -a /cndd3/dburrows/DATA/annotations/rmsk/rmsk.hg38.filt-5ptrim.merge.L1HS_all.canon.sorted.bed \
    > {d}/L1HS_all.bam
    """
    get_ipython().run_cell_magic('bash', '', run)
    # Additional code for each directory can be added here

#this section shows what default arguments will be run if just executing the script
if __name__ == "__main__":
    path = '/cndd/dburrows/DATA/te/rna/PE.counts/ATEM_old/'
    bam_path = '/cndd3/dburrows/DATA/te/rna/PE.bam/'
    dir_list = glob.glob(path + '*Samp*')

    # Create a pool of 10 worker processes
    with multiprocessing.Pool(10) as pool:
        # Map the process_directory function to each item in dir_list
        pool.map(process_directory, dir_list)

    print("All tasks completed.")