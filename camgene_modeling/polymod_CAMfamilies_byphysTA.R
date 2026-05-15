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

# load physiological genotype categorizations
setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/physiology/")
phys <- fromJSON(file = "./physiological_categories_from_TA.json")
setwd("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/")

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
runmodels <- function(genename,tpmdf,degree=5){
  # genename should be a unique gene abbreviation
  # also note, FDR will be run AFTER running ALL models! (i.e. on all genes)
  
  ## get the relevant gene ID
  gid <- cam[match(genename,cam[,'gene_abbr_unique']),'GeneID']
  
  if(gid %in% colnames(tpmdf)){
    # set up base formula
    f <- paste(gid,paste0("~ poly(ZT,",as.character(degree),")"))
    
    # fit the models
    fit0 <- lm(f, data=tpmdf)
    fit1 <- lm(paste0(f,"*phys"), data=tpmdf)
    fit2 <- lm(paste0(f,"*treat"), data=tpmdf)
    fit3 <- lm(paste0(f,"*treat*phys"), data=tpmdf)
    
    adf <- as.data.frame(anova(fit0,fit1,fit2,fit3))
    adf$Model <- c("polyZT","polyZT*treat","polyZT*phys","polyZT*treat*phys")
    adf$GeneID <- gid
    adf$GeneName <- genename
    
    return(adf)
  }
  
}

# run all models; combine results into single df (scaled data)
gnames <- unique(cam$gene_abbr_unique)
dflist <- list()
for(i in 1:length(gnames)){
  print(gnames[i])
  dflist[[i]] <- runmodels(gnames[i],stpm_md,3)
}
allres <- bind_rows(dflist)
allres$FDR_p <- p.adjust(allres$`Pr(>F)`,method = "BH")

# find significant results
sigres <- allres %>% filter(FDR_p<0.05)

# just out of curiosity: try again with degree=5 and see if that makes a difference in the number of significant p-values
dflist <- list()
for(i in 1:length(gnames)){
  print(gnames[i])
  dflist[[i]] <- runmodels(gnames[i],stpm_md,5)
}
d5res <- bind_rows(dflist)
d5res$FDR_p <- p.adjust(d5res$`Pr(>F)`,method = "BH")
sigd5 <- d5res %>% filter(FDR_p<0.05)

# define a function to evaluate the results
format_res <- function(df){
  # add subgenome & gene_family column
  df <- df %>%
    mutate(
      subgenome = sub("_.*", "", GeneName),
      gene_family = GeneName %>%
        sub("^[^_]+_", "", .) %>%
        sub("_[0-9]+$", "", .)
    )
  
}

sigres <- format_res(sigres)
sigd5 <- format_res(sigd5)

# save results for both degrees
write.table(format_res(allres),"./polymod3_scaled_allCAM_byphys_May2026.txt",sep = "\t",col.names = T,row.names = F)
write.table(format_res(d5res),"./polymod5_scaled_allCAM_byphys_May2026.txt",sep = "\t",col.names = T,row.names = F)
