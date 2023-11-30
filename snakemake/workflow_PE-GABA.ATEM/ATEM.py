#CONDA ENV base (python 3.9.12)
#Import packages
#---------------------------------------
import json
import pandas as pd
import sys
import numpy as np
import glob
import pysam

#Import your modules
#---------------------------------------
sys.path.insert(1, '/cndd3/dburrows/CODE/te_rna_ageing/')
import te_rna_f as te
sys.path.insert(1, '/cndd3/dburrows/CODE/admin_tools/')
from admin_tools import admin_functions as adm

#Read in required files for filtering
js = json.load(open(glob.glob('*config*')[0])) #CHANGE TO MAKE MORE FLEXIBLE?
bed_pl = pd.read_csv(js['bed_plus_path'], sep='\t', header=None)
bed_pl.columns =['Chromosome', 'Start', 'End', 'Strand', 'gene_id', 'family_id', 'class_id', 'length', 'full_Start', 'full_End']
bed_pl = bed_pl.drop(columns=['length'])
bed_mi = pd.read_csv(js['bed_minus_path'],sep='\t', header=None)
bed_mi.columns =['Chromosome', 'Start', 'End', 'Strand', 'gene_id', 'family_id', 'class_id', 'length', 'full_Start', 'full_End']
bed_mi = bed_mi.drop(columns=['length'])


# bam_pl = pr.read_bam(snakemake.input.bam_pl, as_df=True) 
# bam_mi = pr.read_bam(snakemake.input.bam_mi, as_df=True) 
bam_pl = te.pysam_subset(snakemake.input.bam_pl) 
bam_mi = te.pysam_subset(snakemake.input.bam_mi) 


#Swap Start + End for minus strand
bam_mi['Start'], bam_mi['End'] = bam_mi['End'], bam_mi['Start']

#File checks
assert sum(bam_pl['Strand'] == '+') == len(bam_pl), 'Some non plus strands assigned to plus bam'
assert sum(bam_mi['Strand'] == '-') == len(bam_mi), 'Some non minus strands assigned to minus bam'
assert sum(bed_pl['Strand'] == '+') == len(bed_pl), 'Some non plus strands assigned to plus bed'
assert sum(bed_mi['Strand'] == '-') == len(bed_mi), 'Some non minus strands assigned to minus bed'

#Define + and - strand files
pl_pars = [bed_pl, bam_pl, snakemake.input.meta_pl, 'plus'] 
mi_pars = [bed_mi, bam_mi, snakemake.input.meta_mi, 'minus'] 
par_list = [pl_pars, mi_pars]


pd.options.mode.chained_assignment = None  # default='warn'
count_df = pd.DataFrame() #empty count matrix
bam_ll = [[],[]] #empty list of lists to store curr_bam indices

# Filter out reads that do not overlap with 5' portion of insertion
for x,par in enumerate(par_list):
    curr_bed = par[0]
    curr_bam = par[1]
    curr_name = pd.read_csv(par[2], sep='\t', header=None)
    assert len(curr_bam) == len(curr_name), 'Bam and metadata files not the same length'
    curr_bam['UMI']=curr_name[0].values #Add UMI column to bam file

    #Loop through each chromosome
    chr_unq = np.unique(curr_bam['Chromosome'].values)
    for i,chr in enumerate(chr_unq):
        print('Aligning to chromosome ' + chr + ' for ' + par[3] + ' strand')
        
        #Slice bed/bam files by chromosome
        chr_bam = curr_bam[curr_bam['Chromosome'] == chr]
        chr_bed = curr_bed[curr_bed['Chromosome'] == chr]
        count_df, bam_ll[x] = te.five_prime_align(chr_bam, chr_bed, count_df, bam_ll[x])

#Add in CPMs as a column
total_reads = pd.read_csv(snakemake.input.n_reads, sep=" ", header=None)[0].values[0] 
count_df['CPM'] = count_df['Count'].values / total_reads * 1000000 

#Save counts matrix and bam indeces
count_df.to_csv(snakemake.output.count_mat, sep='\t', index=False)
np.save(snakemake.output.bam_ind,bam_ll)

#Make txt file of start sites to remove
if len(bam_ll[0]) > 0:
    pl_umi=pd.read_csv(snakemake.input.meta_pl, sep='\t', header=None).iloc[np.setxor1d(np.arange(0,len(bam_pl)) , bam_ll[0].astype(int))]
elif len(bam_ll[0]) == 0:
    pl_umi=pd.read_csv(snakemake.input.meta_pl, sep='\t', header=None)
if len(bam_ll[1]) > 0:
    mi_umi=pd.read_csv(snakemake.input.meta_mi, sep='\t', header=None).iloc[np.setxor1d(np.arange(0,len(bam_mi)) , bam_ll[1].astype(int))]
elif len(bam_ll[1]) == 0:
    mi_umi=pd.read_csv(snakemake.input.meta_mi, sep='\t', header=None)

np.savetxt(snakemake.output.meta_pl_notin,  pl_umi, fmt='%s')
np.savetxt(snakemake.output.meta_mi_notin,  mi_umi, fmt='%s')
