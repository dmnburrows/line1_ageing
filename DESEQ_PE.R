library(DESeq2)
library(tidyverse)

setwd("/cndd/dburrows/DATA/te/rna/PE.counts/DE/")

coldata <- read.csv("ATEM_design.csv", row.names=1)
p_list <- c("ATEM_sub-COUNT.csv", "ATEM_coarse_sub-COUNT.csv")
name_list <- c("granular", "coarse")

for (i in seq_along(p_list)) {
  print(p_list[i])
  print(name_list[i])
  cnts <- read.csv(p_list[i], header=TRUE, check.names=FALSE, row.names=1)
  
  dds <-DESeqDataSetFromMatrix(countData=cnts, 
                               colData=coldata, 
                               design=~AGEYEARS + sex + Cell.Type) 
  
  dds <- DESeq(dds)
  cell_res <- results(dds, alpha=0.05, name="Cell.Type_GLU_vs_GABA")
  write.csv(as.data.frame(cell_res), file=paste(name_list[i],"DESEQ-celltype.csv", sep="_"))

}






