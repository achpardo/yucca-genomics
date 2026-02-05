# Purpose: Run differential expression (DESeq2) between each parent and the respective subgenome of 
## Y. gloriosa.
# Author: Anna Pardo
# Date initiated: Feb. 3, 2026

# import modules
library(readr)
library(DESeq2)
library(dplyr)
library(stringr)

setwd("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/counts/")

# load counts for Y. aloifolia & Y. filamentosa (parental species)
yacounts <- as.data.frame(read_delim("./Ya_counts_withmd.txt",delim = "\t"))
yfcounts <- as.data.frame(read_delim("./YfH1_counts_withmd.txt",delim = "\t"))
# load hybrid's counts
ygcounts <- as.data.frame(read_delim("./Yg_toYgIS_all_correctedmd_countsmatrix.txt",delim = "\t"))

# split ygcounts by subgenome
yfsubg <- ygcounts %>% select(c(sample_name,genotype,time,treat,ZT,species,matches("Yufil")))
yasubg <- ygcounts %>% select(c(sample_name,genotype,time,treat,ZT,species,matches("Yucal")))

# for yacounts & yfcounts: change genotype to = species
yacounts$genotype <- yacounts$species
yfcounts$genotype <- yfcounts$species

# add Ya & Yf data to the relevant subgenomes
allya <- rbind(yacounts,yasubg)
allyf <- rbind(yfcounts,yfsubg)

# save allya & allyf for future sessions
write_delim(allya,file = "./Ya_plusYg-Ya_counts_withmd_2026-02-03.txt",delim = "\t",col_names = T)
write_delim(allyf,file = "./Yf_plusYg-Yf_counts_withmd_2026-02-03.txt",delim = "\t",col_names = T)

# extract a coldata object for each species
cdya <- allya %>% select(c(sample_name,genotype))
cdya$genotype <- factor(cdya$genotype)

cdyf <- allyf %>% select(c(sample_name,genotype))
cdyf$genotype <- factor(cdyf$genotype)

# format count data
rownames(allya) <- allya$sample_name
allya <- allya %>% select(-c(sample_name,genotype,time,treat,ZT,species))
allya <- t(allya)
allya <- as.data.frame(allya)

rownames(allyf) <- allyf$sample_name
allyf <- allyf %>% select(-c(sample_name,genotype,time,treat,ZT,species))
allyf <- as.data.frame(t(allyf))

# run DESeq for Ya vs. Ya subgenome of Yg
dds <- DESeqDataSetFromMatrix(allya,cdya,design = ~ genotype)
deseq <- DESeq(dds)

get_sig_df = function(df,cont){
  df = mutate(df,GeneID = rownames(df),Contrast = cont)
  df = df[which(df$padj<0.05),]
  return(df)
}

# set up reference & comparison levels
refLevels <- "aloifolia"
compLevels <- c()
for (i in unique(cdya$genotype)) {
  if(i[1]!="aloifolia"){
    compLevels <- append(compLevels,i)
  }
}

cont <- list()
for(i in 1:length(compLevels)){
  cont[[i]] <- c("genotype",compLevels[i],refLevels)
}

contrasts <- list()
for(i in 1:length(compLevels)){
  contrasts[[i]] <- paste0(refLevels,"-V-",compLevels[i])
}

dflist <- list()
for(i in 24:length(cont)){
  dflist[[i]] <- get_sig_df(assign(contrasts[[i]],
                          as.data.frame(results(deseq,
                                                contrast=cont[[i]],
                                                alpha=0.05,pAdjustMethod = "fdr"))),contrasts[[i]])
}

resdf <- bind_rows(dflist)
# save results
write_delim(resdf,file = "./DESeqres_Ya_v_YgYasubg.txt",delim = "\t",col_names = T)


# repeat for Yf
dds <- DESeqDataSetFromMatrix(allyf,cdyf,design = ~ genotype)
deseq <- DESeq(dds)

# set up reference & comparison levels
refLevels <- "filamentosa"
compLevels <- c()
for (i in unique(cdyf$genotype)) {
  if(i[1]!="filamentosa"){
    compLevels <- append(compLevels,i)
  }
}

cont <- list()
for(i in 1:length(compLevels)){
  cont[[i]] <- c("genotype",compLevels[i],refLevels)
}

contrasts <- list()
for(i in 1:length(compLevels)){
  contrasts[[i]] <- paste0(refLevels,"-V-",compLevels[i])
}

dflist <- list()
for(i in 1:length(cont)){
  dflist[[i]] <- get_sig_df(assign(contrasts[[i]],
                                   as.data.frame(results(deseq,
                                                         contrast=cont[[i]],
                                                         alpha=0.05,pAdjustMethod = "fdr"))),contrasts[[i]])
}

resdf <- bind_rows(dflist)
# save resdf
write_delim(resdf,file = "./DESeqres_Yf_v_YgYfsubg.txt",delim = "\t",col_names = T)
