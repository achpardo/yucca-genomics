# Purpose: Run polynomial modeling of clock genes by species & physiotype.
# Author: Anna Pardo
# Date initiated: June 24, 2026

library(dplyr)
library(readr)
library(rjson)

# set working directory
setwd("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/")
# load TPM
tpm <- read_delim("./TPM/Yg_toYgIS_allTPM_correctedmd_over1mil.txt",delim = "\t")

# load physiological genotype categorizations
setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/physiology/")
phys <- fromJSON(file = "./physiological_categories_from_TA.json")

# make a 'phys' column
tpm$phys <- as.character(phys[match(tpm$genotype,names(phys))])
# drop rows with NA for phys
tpm <- tpm %>% filter(phys != "NULL")

# load clock gene annotation
clock <- read_csv("../circadian_light_genes_by_orthology_Yucca.csv")
clock <- as.data.frame(clock)

# function to run modeling for a single gene
runmodels <- function(genename,tpmdf,degree=3,annot=clock){
  # genename should be a unique gene abbreviation
  # also note, FDR will be run AFTER running ALL models! (i.e. on all genes)
  
  ## get the relevant gene ID
  gid <- annot[match(genename,annot[,'gene_name_unique']),'GeneID']
  print(gid)
  
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

# run all models; combine results into single df (scaled data)
gnames <- unique(clock$gene_name_unique)
dflist <- list()
for(i in 1:length(gnames)){
  print(gnames[i])
  dflist[[i]] <- runmodels(gnames[i],stpm_md,3)
}
allres <- bind_rows(dflist)
allres$FDR_p <- p.adjust(allres$`Pr(>F)`,method = "BH")

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
write.table(format_res(allres),"./clockgenes_polymod3_scaled_byphys_Jun2026.txt",sep = "\t",col.names = T,row.names = F)
