# Purpose: Regress CAM genes' expression (for each treatment) between each Yg genotype & each parent. 
## This will help determine whether CAM gene expression follows physiological categorizations for Y. gloriosa.
# Author: Anna Pardo
# Date initiated: Mar. 7, 2026

# load modules
library(dplyr)
library(readr)
library(rjson)

# set working directory
setwd("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/")

# load physiological genotype categorizations
phys <- fromJSON(file = "./phys_figures/physiological_CAM_categories.json")

# load TPM
ygtpm <- read_delim("./rna_insilico_genome/TPM/Yg_toYgIS_allTPM_correctedmd_over1mil.txt",delim = "\t")
yatpm <- as.data.frame(read_delim("./rna_reanalysis_30-Jun-2025/TPM/Yalo_TPM_withmd_July2025.txt",delim = "\t"))
yftpm <- as.data.frame(read_delim("./rna_reanalysis_30-Jun-2025/TPM/YfilH1_TPM_withmd_July2025.txt",delim = "\t"))

# for Yg TPM: make a 'phys' column
ygtpm$phys <- as.character(phys[match(ygtpm$genotype,names(phys))])

# load CAM gene annotation
cam <- as.data.frame(read_delim("./rna_insilico_genome/degs_downstream/camgenes_Ya_Yf_orthology_synteny.txt",delim = "\t"))

# subset all three TPM dataframes based on CAM gene IDs only
genecols <- unique(cam$GeneID)
mdcols <- c("sample_name","genotype","treat","ZT","phys")
ygtpm <- ygtpm %>% select(any_of(c(mdcols,genecols)))
ygtpm <- as.data.frame(ygtpm)
yatpm <- yatpm %>% select(any_of(c(mdcols,genecols)))
yftpm <- yftpm %>% select(any_of(c(mdcols,genecols)))

# sort all TPM data by time point & treatment
yatpm <- yatpm %>% arrange(ZT,treat)
yftpm <- yftpm %>% arrange(ZT,treat)
ygtpm <- ygtpm %>% arrange(ZT,treat)

# for ygtpm: filter out certain ZT points
ygtpm <- ygtpm %>% filter(ZT %in% unique(yatpm$ZT))

# define a function to run a regression model for a given gene's expression (all time points for a given treatment)
# another argument: single-genotype TPM for Yg
# calculate the mean of each time point & regress those
single_regression <- function(gid,t,gttpm){
  if(grepl("Yucal", gid)){
    partpm <- yatpm %>% filter(treat==t) %>% select(c(ZT,gid))
    sp <- "Ya"
  } else {
    partpm <- yftpm %>% filter(treat==t) %>% select(c(ZT,gid))
    sp <- "Yf"
  }
  # calculate mean TPM per time point for partpm
  parmeans <- partpm %>% group_by(ZT) %>% summarize_at(vars(gid),list(parental_mean=mean))
  
  # do the same for Yg genotype TPM
  gttpm <- gttpm %>% filter(treat==t)
  geno <- unique(gttpm$genotype)
  gtmeans <- gttpm %>% group_by(ZT) %>% summarize_at(vars(gid),list(genotype_mean=mean))
  
  # merge parmeans & gtmeans such that only the common ZT are kept
  allmeans <- merge(as.data.frame(gtmeans),as.data.frame(parmeans),by="ZT")
  
  model <- lm(allmeans$parental_mean ~ allmeans$genotype_mean)
  info <- list(model,allmeans)
  names(info) <- c("model",paste0(sp,"_",geno,"_means"))
  
  return(info)
}

# test the function
tpm18 <- ygtpm %>% filter(genotype=="18")
testinfo <- single_regression("Yucal.01G165600.v2.1","D",tpm18)

# function to loop through all genotypes & treatments for a single gene
## output dataframe of R2 values, as well as the list of models & means
regress_all <- function(gid){
  alltreat <- unique(ygtpm$treat)
  outlist <- list()
  r2 <- c()
  adjr2 <- c()
  treatment <- c()
  genotype <- c()
  GeneID <- c()
  for(g in 1:length(unique(ygtpm$genotype))){
    outlist[[g]] <- list()
    gttpm <- ygtpm %>% filter(genotype==unique(ygtpm$genotype)[g])
    for(t in 1:length(alltreat)){
      outlist[[g]][[t]] <- single_regression(gid,alltreat[t],gttpm)
      genotype <- append(genotype,unique(ygtpm$genotype)[g])
      treatment <- append(treatment,alltreat[t])
      r2 <- append(r2,summary(outlist[[g]][[t]]$model)$r.squared)
      adjr2 <- append(adjr2,summary(outlist[[g]][[t]]$model)$adj.r.squared)
      GeneID <- append(GeneID,gid)
    }
  }
  
  res <- list(outlist,data.frame(GeneID,genotype,treatment,r2,adjr2))
  names(res) <- c(paste0(gid,"_detailed_results_list"),paste0(gid,"_r2_dataframe"))
  return(res)
}

# loop over all gene IDs to run this
## get out the relevant gene IDs
camgids <- c()
for(c in colnames(ygtpm)){
  if(grepl("Y",c)){
    camgids <- append(camgids,c)
  }
}

all_results <- list()
for(c in 1:length(camgids)){
  if(grepl("Y",camgids[c])){
    all_results[[c]] <- regress_all(camgids[c])
    names(all_results)[[c]] <- camgids[c]
  }
}

# stick all R2 results dataframes together
dflist <- list()
for(i in 1:length(all_results)){
  dflist[[i]] <- as.data.frame(all_results[[i]][[2]])
}
allr2 <- bind_rows(dflist)

# sort allr2 by R2 (highest to lowest)
allr2 <- allr2 %>% arrange(desc(r2))
# add physiology info
allr2$CAMphys <- as.character(phys[match(allr2$genotype,names(phys))])

# add parent of comparison - if gene ID is Yufil parent would be Yf, etc.
par <- c()
for(i in allr2$GeneID){
  if(grepl("Yufil",i)){
    par <- append(par,"Yf")
  } else if(grepl("Yucal",i)) {
    par <- append(par,"Ya")
  }
}
allr2$comp_parent <- par
# add CAM gene families & names
camsub <- cam %>% select(c(GeneID,gene_abbr,gene_abbr_unique))
allr2 <- merge(allr2,camsub,by="GeneID")
allr2 <- allr2 %>% arrange(desc(r2))
# save these results
write_delim(allr2,file = "//wsl$/Ubuntu/home/leviathan22/yucca-genomics/camgene_modeling/Yg_parents_CAMgene_regression_results.txt",
            delim = "\t")
saveRDS(all_results,"//wsl$/Ubuntu/home/leviathan22/yucca-genomics/camgene_modeling/Yg_parents_CAMreg.rds")
