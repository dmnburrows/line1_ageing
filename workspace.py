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

import te_rna_f as ter
sys.path.insert(1, '/cndd3/dburrows/CODE/admin_tools/')
from admin_tools import admin_functions as adm

def multimap(d):
    curr_ = ter.multimap_stats(d+'/L1HS_all.bam')
    curr_.to_csv(d+'/multimap.csv', index=False)
    print('Done ' + d)

#this section shows what default arguments will be run if just executing the script
if __name__ == "__main__":
    path='/cndd/dburrows/DATA/te/rna/PE.counts/l1hs/'
    dir_list = glob.glob(path + '*Samp*')
    
    # Create a pool of 10 worker processes
    with multiprocessing.Pool(10) as pool:
        # Map the process_directory function to each item in dir_list
        pool.map(multimap, dir_list)

    print("All tasks completed.")
    