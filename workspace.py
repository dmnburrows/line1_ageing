#CONDA ENV base_conda (python 3.9.7)
#Import packages
#---------------------------------------
import sys
import os
import glob
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from multiprocessing import Pool
from tqdm import tqdm


#Import your modules
#---------------------------------------
import te_rna_f as ter
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


def my_bedcov(samplepath):
  sample=samplepath.split('/')[-1]
  pref = 'L1'
  bam=f'/cndd3/dburrows/DATA/te/rna/PE.bam/{sample}/Aligned.sortedByCoord.out.bam'
  outdir=f'/cndd/dburrows/DATA/te/rna/PE.genomic.bins/{sample}/'

  sense=f'bedtools coverage -sorted -counts -g /cndd/emukamel/DrachevaLiu_PsychENCODE_SCZ/L1_bins/hg38.genome -bed -s -a /cndd/dburrows/DATA/te/rna/PE.genomic.bins/{pref}_bins.sorted.bed -b {bam} > {outdir}{pref}_bins.{sample}.coverage.sense.bed'
  antisense=f'bedtools coverage -sorted -counts -g /cndd/emukamel/DrachevaLiu_PsychENCODE_SCZ/L1_bins/hg38.genome -bed -S -a /cndd/dburrows/DATA/te/rna/PE.genomic.bins/{pref}_bins.sorted.bed -b {bam} > {outdir}{pref}_bins.{sample}.coverage.antisense.bed'

  if os.path.exists(outdir):
    os.system(sense)
    os.system(antisense)
  else:
    os.mkdir(outdir)
    os.system(sense)
    os.system(antisense)
    
  print(f'Done {sample}')
  return 0


samples=glob.glob('/cndd3/dburrows/DATA/te/rna/PE.bam/Sample*')
#creat pool of parallel worker processes
with Pool(20) as p:
  #imap applies function in parallel
  x=list(tqdm(p.imap(my_bedcov, samples),total=len(samples)))