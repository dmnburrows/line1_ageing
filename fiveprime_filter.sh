#!/bin/bash

# Convert to ensembl format
gtftk convert_ensembl -i rmsk.hg38.subsample.gtf -o rmsk.hg38.subsample.ensemble.gtf

# calculate feature size
gtftk feature_size -i rmsk.hg38.subsample.ensemble.gtf -o rmsk.hg38.subsample.featsize.gtf


# filter by minimum feature size
gtftk select_by_tx_size -i rmsk.hg38.subsample.featsize.gtf -o rmsk.hg38.subsample.filt.gtf -m 600


#Turn into bed file
