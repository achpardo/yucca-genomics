# Purpose: Run DESeq2 on Y. gloriosa RNA-seq data (uniquely mapped to Yg in silico genome).
# Author: Anna Pardo
# Date initiated: Sept. 30, 2025

# import modules
library(readr)
library(DESeq2)
library(dplyr)
library(stringr)

# set working directory (where counts matrix exists)
setwd("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/counts/")

# load counts matrix (with orphan samples, for now)
cm <- read_delim("./Yg_toYgIS_3reps_correctedmd_countsmatrix.txt",delim = "\t")
# relevant function: DESeqDataSetFromMatrix()
# countData: counts matrix with no metadata
# colData: metadata (row names correspond to columns)

# start by setting up countData & colData
# for count data, columns should be samples, rows should be gene IDs
cm <- as.data.frame(cm)
# find columns that are not gene IDs
for(i in colnames(cm)){
  if(startsWith(i,"Yu")==FALSE){
    print(i)
  }
}

countsOnly <- cm %>% select(-c(genotype,time,treat,ZT,species,Total_Reads,condition))
metadata <- cm %>% select(c(sample_name,genotype,time,treat,ZT,species,Total_Reads,condition))

# set sample ID as row names for counts matrix
rownames(countsOnly) <- countsOnly$sample_name
countsOnly <- countsOnly %>% select(-sample_name)

# transpose counts matrix
tcounts <- t(countsOnly)

# try running DESeq with Total_Reads and condition in the formula
dds <- DESeqDataSetFromMatrix(tcounts,metadata,design = ~ Total_Reads + condition)
deseq <- DESeq(dds)

# this takes too much memory to run on local. save files for running on HPC
write.table(tcounts,"./countsmatrix_for_DESeq_withorphans.txt",sep = "\t")
write.table(metadata,"./coldata_forDESeq.txt",sep = "\t",row.names = F)

# test section: try loading in tcounts
test <- read.table("./countsmatrix_for_DESeq_withorphans.txt",sep = "\t",row.names = 1)
md <- read_delim("./coldata_forDESeq.txt",delim = "\t")
