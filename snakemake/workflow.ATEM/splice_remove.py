# Check that spliced reads have been removed correctly
import pysam as sam
import numpy as np
full = sam.AlignmentFile(snakemake.input.splice_bam, 'rb')
nosplice = sam.AlignmentFile(snakemake.input.nosplice_bam, 'rb')

cig_list = []
full_count=0
for x,read in enumerate(full):
    cig_list.append(read.cigarstring)
    full_count+=1
n_splice = len([i for i in cig_list if 'N' in i])
print(str(np.round(len([i for i in cig_list if 'N' in i])/len(cig_list) * 100,decimals=5)) + ' % of reads are spliced')

cig_list = []
filt_count=0
for x,read in enumerate(nosplice):
    cig_list.append(read.cigarstring)
    filt_count+=1
n_splice_filt = len([i for i in cig_list if 'N' in i])

assert n_splice_filt == 0, 'some spliced reads still remain in 5pfilt-tss_nosplice.bam'
assert full_count - n_splice == filt_count, 'not all spliced reads accounted for'
complete='all spliced reads removed'
print(complete)
np.save(snakemake.output.splice_log, complete)
