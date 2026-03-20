# Purpose: Run polynomial models on CAM genes by physiology of genotypes.
# Author: Anna Pardo
# Date initiated: Mar. 19, 2026

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

# load physiological genotype categorizations
phys <- fromJSON(file = "../phys_figures/physiological_CAM_categories.json")

# make a 'phys' column
tpm$phys <- as.character(phys[match(tpm$genotype,names(phys))])
# drop rows with NA for phys
tpm <- tpm %>% filter(phys != "NULL")

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

# set up a function for model running on a given CAM gene
runmodels <- function(genename,tpmdf){
  # genename should be a unique gene abbreviation
  # also note, FDR will be run AFTER running ALL models! (i.e. on all genes)
  
  ## get the relevant gene ID
  gid <- cam[match(genename,cam[,'gene_abbr_unique']),'GeneID']
  
  if(gid %in% colnames(tpmdf)){
    # set up base formula
    f <- paste(gid,"~ poly(ZT,5)")
    
    # fit the models
    fit0 <- lm(f, data=tpmdf)
    fit1 <- lm(paste0(f,"*phys"), data=tpmdf)
    fit2 <- lm(paste0(f,"*treat"), data=tpmdf)
    fit3 <- lm(paste0(f,"*treat*phys"), data=tpmdf)
    
    adf <- as.data.frame(anova(fit0,fit1,fit2,fit3))
    adf$Model <- c("polyZT","polyZT*phys","polyZT*treat","polyZT*treat*phys")
    adf$GeneID <- gid
    adf$GeneName <- genename
    
    return(adf)
  }
  
}

# run all models; combine results into single df (unscaled data)
dflist <- list()
gnames <- unique(cam$gene_abbr_unique)
for(i in 1:length(gnames)){
  print(gnames[i])
  dflist[[i]] <- runmodels(gnames[i],tpm)
}
allres_unscaledtpm <- bind_rows(dflist)
allres_unscaledtpm$TPM_Type <- "raw"

# repeat for scaled data
dflist <- list()
for(i in 1:length(gnames)){
  print(gnames[i])
  dflist[[i]] <- runmodels(gnames[i],stpm_md)
}
allres_scaled <- bind_rows(dflist)
allres_scaled$TPM_Type <- "scaled"

# stick these two results together
allres <- rbind(allres_scaled,allres_unscaledtpm)

# calculate FDR
allres$FDR_p <- p.adjust(allres$`Pr(>F)`,method = "BH")

# find significant results
sigres <- allres %>% filter(FDR_p<0.05)

# actually: results do not change between scaled & unscaled TPM. run FDR independently on each of these and save both. I will pick one to use in future
allres_scaled$FDR_p <- p.adjust(allres_scaled$`Pr(>F)`,method = "BH")
allres_unscaledtpm$FDR_p <- p.adjust(allres_scaled$`Pr(>F)`,method = "BH")
write.table(allres_scaled,"./polymodres_scaled_allCAM_byphys.txt",sep = "\t",col.names = T,row.names = F)
write.table(allres_unscaledtpm,"./polymodres_rawTPM_allCAM_byphys.txt",sep = "\t",col.names = T,row.names = F)
