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


def process_directory(s):
    run = f""" 
    samtools view {bam_path}/{s}/Aligned.sortedByCoord.out.bam | wc -l > {path}/{s}/total_reads.txt
    samtools view -h -F 0x10 {bam_path}/{s}/Aligned.sortedByCoord.out.bam > {path}/{s}/plus.bam
    samtools view -h -f 0x10 {bam_path}/{s}/Aligned.sortedByCoord.out.bam > {path}/{s}/minus.bam
    samtools view -b -h -L /cndd3/dburrows/DATA/annotations/rmsk/rmsk.hg38.filt-5ptrim.plus.full.sort.bed {path}/{s}/plus.bam > {path}/{s}/plus.5pfilt_withsplice.bam
    samtools index {path}/{s}/plus.5pfilt_withsplice.bam
    samtools view -b -h -L /cndd3/dburrows/DATA/annotations/rmsk/rmsk.hg38.filt-5ptrim.minus.full.sort.bed {path}/{s}/minus.bam > {path}/{s}/minus.5pfilt_withsplice.bam
    samtools index {path}/{s}/minus.5pfilt_withsplice.bam
    
    """
    get_ipython().run_cell_magic('bash', '', run)
    

    print('Done' + s)
    

#this section shows what default arguments will be run if just executing the script
if __name__ == "__main__":
    bam_path = '/cndd3/dburrows/DATA/te/rna/PE.bam/'
    path = '/cndd/dburrows/DATA/te/rna/PE.counts/ATEM/'
    samp_list = ['Sample_1134_GABA','Sample_1539_GABA', 'Sample_1648_GABA', 'Sample_198034-1435_GABA', 'Sample_4369_GABA',
    'Sample_4414_GABA', 'Sample_4545_GABA', 'Sample_5077_GABA', 'Sample_5387_GABA', 'Sample_5401_GABA',
    'Sample_5446_GABA', 'Sample_6007_GABA', 'Sample_HCT15HBNA032_GABA', 'Sample_HCT16HCQA020_GABA', 
    'Sample_HCT16HECA028_GABA']
    
    # Create a pool of 10 worker processes
    with multiprocessing.Pool(10) as pool:
        # Map the process_directory function to each item in dir_list
        pool.map(process_directory, samp_list)

    print("All tasks completed.")
    
