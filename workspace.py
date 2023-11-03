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


def my_bedcov(samplepath, qthresh=30):
  sample=samplepath.split('/')[-1]
  pref = 'L1'
  bam=f'/cndd3/dburrows/DATA/te/rna/PE.bam/{sample}/Aligned.sortedByCoord.out.bam'
  outdir=f'/cndd/dburrows/DATA/te/rna/PE.genomic.bins/{sample}/'
  
  if os.path.exists(outdir) ==False:
    os.mkdir(outdir)
    
  outfile=f'{pref}_bins.{sample}.q{qthresh}.sense.coverage.bed'
  #filter multimapped reads, and also get reads pair aligning to sense strand  
  cmd=f'/usr/bin/samtools view -q {qthresh} -M -f64 -b -L /cndd/dburrows/DATA/te/rna/PE.genomic.bins/{pref}_bins.merged.bed {bam} | '
  cmd+=f' bedtools coverage -sorted -c -s -a /cndd/dburrows/DATA/te/rna/PE.genomic.bins/{pref}_bins.sorted.bed -b - > {outdir}{outfile} '
  os.system(cmd)
    
  cmd=f'/usr/bin/samtools view -q {qthresh} -M -f128 -b -L /cndd/dburrows/DATA/te/rna/PE.genomic.bins/{pref}_bins.merged.bed {bam} | '
  cmd+=f' bedtools coverage -sorted -c -S -a /cndd/dburrows/DATA/te/rna/PE.genomic.bins/{pref}_bins.sorted.bed -b - >> {outdir}{outfile} '
  os.system(cmd)

  outfile=f'{pref}_bins.{sample}.q{qthresh}.antisense.coverage.bed'
  cmd=f'/usr/bin/samtools view -q {qthresh} -M -f64 -b -L /cndd/dburrows/DATA/te/rna/PE.genomic.bins/{pref}_bins.merged.bed {bam} {bam} | '
  cmd+=f' bedtools coverage -sorted -c -S -a  /cndd/dburrows/DATA/te/rna/PE.genomic.bins/{pref}_bins.sorted.bed -b - > {outdir}{outfile} '
  os.system(cmd)
  cmd=f'/usr/bin/samtools view -q {qthresh} -M -f128 -b -L /cndd/dburrows/DATA/te/rna/PE.genomic.bins/{pref}_bins.merged.bed {bam} | '
  cmd+=f' bedtools coverage -sorted -c -s -a /cndd/dburrows/DATA/te/rna/PE.genomic.bins/{pref}_bins.sorted.bed -b - >> {outdir}{outfile} '
  os.system(cmd)
    
  print(f'Done {sample}')
  return 0


samples=glob.glob('/cndd3/dburrows/DATA/te/rna/PE.bam/Sample*')
#creat pool of parallel worker processes
with Pool(20) as p:
  #imap applies function in parallel
  x=list(tqdm(p.imap(my_bedcov, samples),total=len(samples)))