import numpy as np
import pandas as pd
import pyranges as pr

def check(check_bool, check_bam, sign):
    bam_ = np.char.add((check_bam[check_bam['Strand'] == sign]['Start'].values).astype(str),  (check_bam[check_bam['Strand'] == sign]['Chromosome'].values).astype(str))
    bool_ = np.char.add(check_bool['Start'].values.astype(str), check_bool['Chromosome'].values.astype(str) )
    assert sum(np.in1d(bam_, bool_)) == len(bam_), 'Merged BAM file is missing some strand reads'
    assert sum(np.in1d(bool_, bam_)) == len(bool_), 'Merged BAM file is missing some strand reads'


#Final assertion
out_bam = pr.read_bam(snakemake.input.merged_bam, as_df=True) 
bam_ll = np.load(snakemake.input.bam_ind, allow_pickle=True)
bam_pl = pr.read_bam(snakemake.input.bam_pl, as_df=True) 
bam_mi = pr.read_bam(snakemake.input.bam_mi, as_df=True) 

assert len(out_bam[out_bam['Strand'] == '+']) == len(bam_ll[0]), 'Merged BAM file is missing some plus strand reads'
assert len(out_bam[out_bam['Strand'] == '-']) == len(bam_ll[1]), 'Merged BAM file is missing some minus strand reads'
check(bam_pl.iloc[bam_ll[0]],out_bam, '+')
check(bam_mi.iloc[bam_ll[1]],out_bam, '-')
complete='All reads accounted for'
print(complete)
np.save(snakemake.output.log_complete, complete)
