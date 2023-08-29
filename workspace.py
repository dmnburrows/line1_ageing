
#CONDA ENV base (python 3.9.7)
#Import packages
#---------------------------------------
import sys
import os
import glob
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from Bio import SeqIO
import pyranges as pr
import scanpy as sc
import pysam as sam

#Import your modules
#---------------------------------------
import te_rna_f as te
sys.path.insert(1, '/cndd3/dburrows/CODE/admin_tools/')
from admin_tools import admin_functions as adm

# Define paths
#----------------------------------------------------------------------
l_code = '/Users/dominicburrows/Dropbox/PhD/Analysis/my_scripts/GitHub/'
l_data = '/Users/dominicburrows/Dropbox/PhD/analysis/Project/'
l_fig = '/Users/dominicburrows/Dropbox/PhD/figures/'

s_code = '/cndd3/dburrows/CODE/'
s_data = '/cndd3/dburrows/DATA/'
s_fig = '/cndd3/dburrows/FIGS/'

sys.version

curr_l = glob.glob('/cndd/dburrows/DATA/te/rna/PE.counts/ATEM/' + 'Samp*')

for x,c in enumerate(curr_l):
    df = te.multimap_stats(c+'/5pfilt-tss.bam')
    #save csv 
    df.to_csv(c + '/multimap_stats.csv')
    print(str(x) + ' of ' + str(len(curr_l)))
    
