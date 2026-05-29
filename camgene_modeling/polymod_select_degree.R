# Purpose: Run polynomial models on CAM genes by physiology of genotypes (physiotypes from TA data).
# Author: Anna Pardo
# Date initiated: May 14, 2026

library(dplyr)
library(readr)
library(rjson)
# set working directory
setwd("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/")
# load TPM
tpm <- read_delim("./TPM/Yg_toYgIS_allTPM_correctedmd_over1mil.txt",delim = "\t")

for(c in colnames(tpm)){
  if(startsWith(c,"Yu")==FALSE){
    print(c)
  }
}

# scale TPM before running models
stpm <- as.data.frame(tpm)
row.names(stpm) <- stpm$sample_name
stpm <- stpm %>% select(-c(sample_name,genotype,time,treat,ZT,species,phys))
stpm <- scale(stpm)
stpm <- as.data.frame(stpm)
# add back metadata
stpm$sample_name <- row.names(stpm)
md <- as.data.frame(tpm %>% select(c("sample_name","genotype","phys","treat","ZT")))
stpm_md <- merge(md,stpm)

# load CAM gene annotation
cam <- as.data.frame(read_delim("./degs_downstream/camgenes_Ya_Yf_orthology_synteny.txt",delim = "\t"))
# load lists of time-structured genes from Yg (to select some for analysis that aren't CAM)
setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/masigpro/")
camts <- as.data.frame(read_delim("./CAM_TSgenes_9gt.txt",delim = "\t"))
c3ts <- as.data.frame(read_delim("./C3+CAM_TSgenes_9gt.txt",delim = "\t"))
facts <- as.data.frame(read_delim("./facCAM_TSgenes_9gt.txt",delim = "\t"))
tsgenes <- bind_rows(camts,c3ts,facts)
tsgenes <- tsgenes %>% filter(!GeneID %in% unique(cam$GeneID))
tsgenes <- unique(tsgenes$GeneID)

# select 10 random time-structured genes (no replacement) that are not CAM genes
testgenes <- sample(tsgenes,10)

# set up a function to run models of varying degrees on a given gene
runmod <- function(gid,tpmdf=stpm_md){
  if(gid %in% colnames(tpmdf)){
    # fit the models
    fitx <- lm(paste(gid,"~ poly(ZT,1)"), data=tpmdf)
    fit0 <- lm(paste(gid,"~ poly(ZT,2)"), data=tpmdf)
    fit1 <- lm(paste(gid,"~ poly(ZT,3)"), data=tpmdf)
    fit2 <- lm(paste(gid,"~ poly(ZT,4)"), data=tpmdf)
    fit3 <- lm(paste(gid,"~ poly(ZT,5)"), data=tpmdf)
    fit4 <- lm(paste(gid,"~ poly(ZT,6)"), data=tpmdf)
    
    adf <- as.data.frame(anova(fitx,fit0,fit1,fit2,fit3,fit4))
    adf$Model <- c("ZT1","ZT2","ZT3","ZT4","ZT5","ZT6")
    adf$GeneID <- gid
    
    return(adf)
  } else {
    return("Error: gene not in TPM!")
  }
  
}

dflist <- list()
for(i in 1:length(testgenes)){
  dflist[[i]] <- runmod(testgenes[i])
}
allres <- bind_rows(dflist)
# calculate FDR
allres$FDR_p <- p.adjust(allres$`Pr(>F)`,method = "BH")
allres <- allres %>% filter(Model != "ZT1")
sigres <- allres %>% filter(FDR_p < 0.05)
sigres %>% group_by(Model) %>% summarize(count=n()) # count number of genes significant for each model

# find which model has the lowest p-value for each gene
for (i in unique(allres$GeneID)) {
  df <- allres %>% filter(GeneID==i)
  df <- df %>% arrange(FDR_p)
  print(paste0("Model with lowest p for ",i," : ",df$Model[1]))
}

