# Purpose: Extract results from ASE re-run.
# Author: Anna Pardo
# Date initiated: Apr. 10, 2026

# import modules
library(readr)
library(DESeq2)
library(dplyr)
library(stringr)

# set working directory
setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/differential_expression")
# load RData generated on the HPC
load("./ase_gttreat_pre_resextract.RData")

# extract results
df_list <- list()
for(i in 1:length(compare_list)){
  print(compare_list[i][[1]])
  condition_list <- unlist(strsplit(as.character(compare_list[i][[1]]),"-V-"))
  contrast <- paste0(condition_list[1],"-V-",condition_list[2])
  print(contrast)
  df_list[[i]] = assign(paste0(contrast,"_sig"),get_sig_df(assign(contrast,as.data.frame(results(dds,contrast=c("gtsg",condition_list[2],condition_list[1]),alpha=0.05,pAdjustMethod = "fdr"))),contrast))
}
DE_Gene_df = bind_rows(df_list)
write_delim(DE_Gene_df,outfile,delim="\t")