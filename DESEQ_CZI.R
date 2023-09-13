library(DESeq2)
library(tidyverse)

setwd("/cndd3/dburrows/DATA/te/rna/CZI.counts/DESEQ/")
coldata <- read.csv("ATEM_design.csv", row.names=1)

#Loop over all the files in the current directory
files <- list.files(pattern="*.ATEM.csv")

#Loop over each csv in files and run DESEQ
for(i in files) {                                          
  input <- i      
    cnts <- read.csv(input, header=TRUE, check.names=FALSE, row.names=1)
    dds <-DESeqDataSetFromMatrix(countData=cnts, 
                                colData=coldata, 
                                design=~Subject.Age+Subject.Sex) 
    dds <- DESeq(dds)
    res <- results(dds, alpha=0.1)
    output <- paste0(substr(i, 1, nchar(i)-8), ".DESEQ_sex.csv")
    print(output)
    write.csv(as.data.frame(res), file=output)

}









