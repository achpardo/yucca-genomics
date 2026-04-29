# Purpose: Use DESeq2 to determine homeolog expression bias (HEB) in Y. gloriosa.
# Author: Anna Pardo
# Date initiated: April 6, 2026

# import modules
library(readr)
library(DESeq2)
library(dplyr)
library(stringr)

# set working directory
setwd("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/counts/")

# load counts & colData
counts <- read_delim("./Yg_toYgIS_syntelog_counts_forHEB_DESeq.txt",delim = "\t")
coldata <- read_delim("./coldata_deseq_heb.txt",delim = "\t")

# set syntelog IDs as row IDs
counts <- as.data.frame(counts)
row.names(counts) <- counts$syntelogID
counts <- counts %>% select(-syntelogID)

# in coldata: set up genotype_treatment column
coldata$gttreat <- paste(coldata$genotype,coldata$treat,sep = "_")
coldata$gtsg <- paste(coldata$gttreat,coldata$subgenome,sep = "_")

# find number of reps for each of these conditions (column: gtsg)
repcounts <- coldata %>% group_by(gtsg) %>% count()
repcounts %>% filter(n<3)
# there are no conditions with <3 replicates - good.

# drop some information from coldata
cd <- coldata %>% select(c(sample_name,gtsg))

# run DESeq by gtsg
dds <- DESeqDataSetFromMatrix(counts,cd,design = ~ gtsg)
dds <- DESeq(dds)

# set up reference levels (desired: Ya as reference for all genotype_treatment combinations)
refLevels <- c()
for(i in unique(coldata$gtsg)){
  if(grepl("aloifolia",i)){
    refLevels <- append(refLevels,i)
  }
}
