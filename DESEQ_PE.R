library(DESeq2)
library(tidyverse)

setwd("/cndd/dburrows/DATA/te/rna/PE.counts/DESEQ/")

coldata <- read.csv("ATEM_design.csv", row.names=1)
cnts <- read.csv("ATEM_sub-COUNT.csv", header=TRUE, check.names=FALSE, row.names=1)

dds <-DESeqDataSetFromMatrix(countData=cnts, 
                             colData=coldata, 
                             design=~Cell.Type+sex+race+AGEYEARS) 

dds <- DESeq(dds)
age_res <- results(dds, alpha=0.1)
cell_res <- results(dds, alpha=0.1, name="Cell.Type_GLU_vs_GABA")
sex_res <- results(dds, alpha=0.1, name="sex_Male_vs_Female")
race_res <- results(dds, alpha=0.1, name="race_White_vs_Black")

write.csv(as.data.frame(age_res), file="DESEQ_age.csv")
write.csv(as.data.frame(cell_res), file="DESEQ_cell.csv")
write.csv(as.data.frame(sex_res), file="DESEQ_sex.csv")
write.csv(as.data.frame(race_res), file="DESEQ_race.csv")



