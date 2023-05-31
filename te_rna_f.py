import numpy as np

#===================================================
def find_intersect(bam_remaining, flat_bed_pos, flat_bed_ind):
#===================================================

    """
    This function finds the intersection between a bam file and a flattened bed file of insertion positions, 
    and returns the UMIs, bed indeces, and bam indeces of the reads that overlap.

    Inputs:
    bam_remaining: bam file of reads that have not yet been counted
    chr_bed: bed file of insertions on a given chromosome
    flat_bed_pos: flattened vector of 5' insertion positions across all insertions in chr_bed
    flat_bed_ind: flattened vector of indeces for each region that maps it back onto the original chr bed file

    Outputs:
    umi: vector of UMIs that overlap with flattened bed
    bedind: vector of indeces in original bed file where umi_v reads have aligned
    ind: vector of pd row indeces of original bam file where reads have aligned

    """
    

    _int = np.intersect1d(bam_remaining['Start'].values, flat_bed_pos, return_indices=True)  #Find indeces (in the bam file of 5' aligned reads only) of reads whose tss overlaps with flattened bed vector
    umi = bam_remaining['UMI'].iloc[_int[1]].values #vector of UMIs that overlap with flattened bed
    bedind = flat_bed_ind[_int[2]] #vector of indeces in original bed file where umi_v reads have aligned
    ind = bam_remaining.index[_int[1]].values #vector of indeces of original bam file where reads have aligned
    assert len(umi) == len(bedind), 'Bam and bed slices not the same length'

    return(umi, bedind, ind)



#=========================================
def five_prime_align(chr_bam, chr_bed,  count_df, bam_ind):
#=========================================
    """
    This function takes a bam file of reads and a bed file of TE insertions and
    filters out reads that do not overlap with the 5' portion of the TE insertion. 
    The output is a csv file of counts at each TE insertion. Multi-read UMIs have
    normalised counts by dividing by the number of multi-reads for each UMI.

    Inputs:
        chr_bam: bam file of reads aligned to a single chromosome
        chr_bed: bed file of TE insertions on a single chromosome
        count_df: dataframe of counts at each TE insertion
        bam_ind: vector of indeces from original full bam file of reads that have aligned to 5' ends at their tss

    Outputs:
        count_df: dataframe of counts at each TE insertion
        bam_ind: vector of indeces from original full bam file of reads that have aligned to 5' ends at their tss
    """

    import numpy as np
    import pandas as pd
    import os
    import sys
    import te_rna_f as te
    sys.path.insert(1, '/cndd3/dburrows/CODE/admin_tools/')
    from admin_tools import admin_functions as adm

    # Filter out reads that do not overlap with 5' portion of insertion

    #Generate flattened vector of 5' insertion positions and their indeces
    all_bed_pos = np.asarray([(np.arange(chr_bed['Start'].values[i], chr_bed['End'].values[i]+1), np.full(chr_bed['End'].values[i]+1 - chr_bed['Start'].values[i],i)) for i in range(len(chr_bed))]) #list of all bed positions for each insertion and their indeces
    assert len(all_bed_pos) == len(chr_bed), 'Not all bed positions accounted for'

    flat_bed_pos = np.ravel(np.asarray(all_bed_pos)[:,0,:]) # flattened vector of all 5' regions across all insertions
    flat_bed_ind = np.ravel(np.asarray(all_bed_pos)[:,1,:]) # flattened vector of indeces for each region that maps it back onto the original bed file
    assert len(flat_bed_pos) == len(flat_bed_ind), 'Bed position and index vectors not the same length'

    # get BAM file of final aligning reads
    bam_bool = np.in1d(chr_bam['Start'].values, flat_bed_pos) #Boolean of indeces of reads whose tss overlaps with bed files
    bam_final = chr_bam[bam_bool] #Final bam reads that have aligned to 5' ends


    bam_remaining = bam_final

    #Mop up reads and their locations until none remaining
    umi_v, bedind_v, ind_v = [],[],[]
    #Loop until all reads have been accounted for
    while len(bam_remaining) > 0:
        umi, bedind, ind = te.find_intersect(bam_remaining, flat_bed_pos, flat_bed_ind)
        umi_v = np.append(umi_v, umi) #UMIs of reads that have aligned to 5' ends
        bedind_v=np.append(bedind_v,bedind) #Indeces of chr_bed insertions where each UMI has aligned
        ind_v = np.append(ind_v, ind) # pandas row indeces of chr_bam file where each UMI comes from
        bam_remaining = bam_remaining.drop(ind) #Drop reads that have already been counted
    if len(bedind_v) > 0: bedind_v = bedind_v.astype(int)
    assert len(bam_final) == len(umi_v) == len(bedind_v) == len(ind_v), 'Not all reads accounted for'
    assert len(np.unique(ind_v)) == len(ind_v), 'Some reads counted twice'
    #append index vector for creation of bam file of reads that have aligned to 5' ends
    if len(bedind_v) > 0: bam_ind = np.append(bam_ind, ind_v.astype(int))

    #Sort by UMI
    sort_umi, sort_bedind = adm.sort_2list(umi_v, bedind_v)

    #Get counts!
    #==============================================================================
    slice_umi, slice_count= np.unique(sort_umi, return_counts=True) #Find unique multi-UMIs and their counts
    counts_v = np.ones(len(sort_umi))

    #Loop over multi-UMIs
    for s in slice_umi[np.where(slice_count > 1)]:
        sort_umi_ind = np.ravel(np.argwhere(np.in1d(sort_umi,s))) #Indeces of sorted UMI list where current multi-UMI is found

        #Assert
        #Check that all identified multi-UMIs have at least 2 reads
        if len(sort_umi_ind) < 2:
            assert False, 'Identified multi-UMIs have only 1 or 0 reads'

        assert sum(counts_v[sort_umi_ind] != 1) == 0, str(x) #'This multi-UMIs has already been counted'
        counts_v[sort_umi_ind] = 1/len(sort_umi_ind) 

    assert len(sort_umi) ==  len(sort_bedind) == counts_v.shape[0], 'Output vectors wrong shape'

    curr_df = chr_bed.iloc[sort_bedind]
    curr_df['Count'] = counts_v 
    count_df = pd.concat([count_df, curr_df]) #Add them to the count matrix
    return(count_df, bam_ind)




#==============================================================================
def calculate_age(milli_div, subsitution_rate=2.2):
#==============================================================================
    
    p = milli_div / 1000  # The milliDiv column in the `rmsk.txt` file.
    p_part = (4 / 3) * p
    jc_dist = -0.75 * (np.log(1 - p_part))
    (jc_dist * 100) / (subsitution_rate * 2 * 100) * 1000
    return jc_dist

#=========================================
def read_rmsk(filename: str):
#=========================================

    import gzip
    from math import log
    import pandas as pd

    # read first line to check if it is a valid rmsk file
    if filename.endswith(".gz"):
        with gzip.open(filename) as f:
            line = f.readline()
    else:
        with open(filename) as f:
            line = f.readline()

        assert (
            line
            == "   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat\n"
        ), "Not a valid rmsk file"

    # setup converter functions
    strand_conv = lambda x: "-" if x == "C" else "+"
    coord_conv = lambda x: int(x.rstrip(")").lstrip("("))
    perc_conv = lambda x: float(x) * 10

    convs = {
        "milliDiv": perc_conv,
        "milliDel": perc_conv,
        "milliIns": perc_conv,
        "genoLeft": coord_conv,
        "strand": strand_conv,
        "repStart": coord_conv,
        "repLeft": coord_conv,
    }

    # read the rmsk file
    df = pd.read_csv(
        filename,
        skiprows=3,
        delim_whitespace=True,
        names=[
            "swScore",
            "milliDiv",
            "milliDel",
            "milliIns",
            "genoName",
            "genoStart",
            "genoEnd",
            "genoLeft",
            "strand",
            "repName",
            "repClassFamily",
            "repStart",
            "repEnd",
            "repLeft",
            "id",
        ],
        converters=convs,
    )

    # split repClassFamily into repClass and repFamily on /
    df[["repClass", "repFamily"]] = df["repClassFamily"].str.split("/", expand=True)
    df.drop("repClassFamily", axis=1, inplace=True)

    # calculate length of each repeat
    df["length"] = df.apply(
        lambda x: x["repEnd"] - x["repLeft"]
        if x["strand"] == "-"
        else x["repEnd"] - x["repStart"],
        axis=1,
    )

    # calculate age of each repeat
    df["age"] = df["milliDiv"].apply(calculate_age)

    return df