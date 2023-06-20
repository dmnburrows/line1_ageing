#!/bin/bash

snakemake --use-conda --cores 16 --configfile config.PE-GLU.json >& /cndd3/dburrows/DATA/te/rna/PE.counts/ATEM/log.snakemake_run