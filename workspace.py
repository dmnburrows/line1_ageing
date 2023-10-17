#CONDA ENV base_conda (python 3.9.7)
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

%load_ext autoreload
sys.version


df = pd.read_csv('/cndd/dburrows/DATA/te/rna/PE.counts/DE/ATEM_CPM.csv', index_col=0)
meta = pd.read_csv('/datasets/Public_Datasets/Dracheva_PsychEncode_development/processed/metadata_RNA_QC.tsv.gz', sep='\t') 
meta = meta[meta['RNA_passQC']]

#Filter low CPM data
thresh = 1 
ind = np.mean(df, axis=1) >thresh #indeces of ones to keep
sub_df = df.loc[ind]

#load in metadata
meta_ = meta
INF = meta_[meta_['age'] < 2]# <2 infancy
ECH = meta_[meta_['period'] == 'Early_Childhood'] # 2-5 early childhood
ECH = ECH[ECH['age'] >= 2]
LCH = meta_[meta_['period'] == 'Late_Childhood'] #5-12 late childhood
ADO = meta_[meta_['period'] == 'Adolescence'] #12-20 adolescence
ADU = meta_[meta_['period'] == 'Adulthood'] #20-50 adulthood
LADU = meta_[meta_['period'] == 'Late_Adulthood'] #50-80 late adulthood

#Generate dictionary of null distributions for all comparisons -> U statistics only
#=============
gene_mat = pd.read_csv('/datasets/Public_Datasets/Dracheva_PsychEncode_development/processed/rna_seq/allSamples_rsem_genes_results_TPM_mod.txt.gz', sep='\t', header=0, index_col=0)
#filter low tpm data
thresh = 1
ind = np.mean(gene_mat, axis=1) >thresh #indeces of ones to keep
sub_df = gene_mat.loc[ind]

#make into df with RNA, age, class, celltype
meta_ = meta[['celltype', 'sex', 'age', 'period']]
meta_['period'][meta_['age'] < 2] = 'Infancy'
curr_df = sub_df[meta['sample'].values]

if mode == 'coarse':
    group_df = {'period':[], 'celltype':[], 'Class':[], 'RNA':[]}
    
    #make into df
    group_df['RNA'] = np.ravel(curr_df)
    group_df['Class'] = np.repeat(curr_df.index.values,len(curr_df.columns.values))
    group_df['period'] = meta_['period'].tolist()*len(curr_df.index)
    group_df['celltype'] = [i.split('_')[1] for i in curr_df.columns.values.tolist()*len(curr_df.index)]
    group_df = pd.DataFrame(group_df)

class_l = group_df['Class'].unique()
cell_l = 'GLU', 'GABA'
period_l = 'Early_Childhood', 'Late_Childhood', 'Adolescence', 'Adulthood', 'Late_Adulthood'
comp_df = ter.inf_paired_comp(class_l, cell_l, period_l, group_df)
comp_df.to_csv('/cndd/dburrows/DATA/te/rna/PE.counts/DE/ATEM_null_test.csv', header=True, index=False)