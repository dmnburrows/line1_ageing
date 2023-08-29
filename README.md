# te_ageing

This repo is for the preprocessing and analysis transposable elements across the lifepsan in aligned genomic data (WGS, rna, methylation)

## What is this repo for?
* the processing of RNA datasets from human brain donors
* the construction of algorithms for quantifying active TE RNA expression
* the analysis of transposable element RNA expression across the lifespan


## What does this repo contain?
* Modules contain functions for TE analysis
* Accompanying ipynotebooks demonstrate how to use the modules
* shell scripts for running bioinformatics pipelines
* directories for running snakemake with ATEM


### Modules
'te_rna_f.py' - code for analysing TE RNA expression and for Active Transposable Element Mapping (ATEM)

'splitbam_by_cluster.py' - splitting bam file into different clusters

'DESEQ_PE.R' - R script for running DESEQ on bulk GLU + GABA datasets

'DESEQ_CZI.R' - R script for running DESEQ on psueodbulk single cell neuronal data


### Notebooks

'te_rna_filtering.ipynb' - different strategies for filtering RNA reads to find active TEs - implementation of Active Transposable Element Mapping (ATEM)

'te_rna_analysis.ipynb' - analysis of ATEM across lifespan and cell types
