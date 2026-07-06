# Purpose: Run polynomial modeling on genes from various pathways of interest (photorespiration, 
## nitrogen metabolism, and citrate cycle).
# Author: Anna Pardo
# Date initiated: June 25, 2026

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
# drop columns with all NA values
stpm_md <- stpm_md %>% 
  select(
    where(
      ~sum(!is.na(.x)) > 0
    )
  )


# load pathway annotation
pathgenes <- as.data.frame(read_csv(file = "../photresp_N_citrate_genes_Yucca.csv"))

# function to run modeling for a single gene
runmodels <- function(genename,tpmdf,degree=3,annot=clock){
  # genename should be a unique gene abbreviation
  # also note, FDR will be run AFTER running ALL models! (i.e. on all genes)
  
  ## get the relevant gene ID
  if("gene_name_unique" %in% colnames(annot)){
    gid <- annot[match(genename,annot[,'gene_name_unique']),'GeneID']
  } else {
    gid <- annot[match(genename,annot[,'gene_abbr_unique']),'GeneID']
  }
  
  if(gid %in% colnames(tpmdf)){
    print(gid)
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

gnames <- unique(pathgenes$gene_name_unique)
dflist <- list()
for(i in 1:length(gnames)){
  print(gnames[i])
  dflist[[i]] <- runmodels(gnames[i],stpm_md,3,pathgenes)
}
allres <- bind_rows(dflist)
allres$FDR_p <- p.adjust(allres$`Pr(>F)`,method = "BH")

# repeat for PPDK-RPs
## load new CAM annotation
setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/")
cam <- read_delim("./camgenes_Ya_Yf_orthology_synteny.txt",delim = "\t")
rp <- as.data.frame(cam %>% filter(gene_abbr == "PPDK-RP"))
rpgenes <- unique(rp$gene_abbr_unique)

dflist <- list()
for(i in 1:length(rpgenes)){
  print(rpgenes[i])
  dflist[[i]] <- runmodels(rpgenes[i],stpm_md,3,rp)
}
rpres <- bind_rows(dflist)
rpres$FDR_p <- p.adjust(rpres$`Pr(>F)`,method = "BH")

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

allres <- format_res(allres)
rpres <- format_res(rpres)

# load results for other gene sets (to run FDR on everything together)
clockres <- read_delim("./physiology/clockgenes_polymod3_scaled_byphys_Jun2026.txt",delim = "\t")
camres <- read_delim("./masigpro/polymod3_scaled_allCAM_byphys_May2026.txt",delim = "\t")

# put all data together
datalist <- list(allres,rpres,clockres,camres)
totalres <- bind_rows(datalist)
# re-run FDR
totalres$FDR_p_largen <- p.adjust(totalres$`Pr(>F)`,method = "BH")
# save results
write.table(totalres,file = "./polymod3_res_all_pathways_Jun2026.txt",sep = "\t",col.names = T,row.names = F)
