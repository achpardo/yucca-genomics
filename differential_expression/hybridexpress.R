# Purpose: Try out HybridExpress for gene categorization.
# Author: Anna Pardo
# Date initiated: Feb. 11, 2026

library(HybridExpress)
library(DESeq2)
library(readr)
library(dplyr)
library(SummarizedExperiment)

setwd("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/")

##### Set up SummarizedExperiment object #####

# load counts matrix for each species (with syntelog information, if already present)
ygcounts <- read_delim("./counts/Yg_toYgIS_syntelog_counts_forHEB_DESeq.txt",delim = "\t")
ygcounts <- as.data.frame(ygcounts)

yacounts <- as.data.frame(read_delim("./counts/Ya_counts_withmd.txt",delim = "\t"))
yfcounts <- as.data.frame(read_delim("./counts/YfH1_counts_withmd.txt",delim = "\t"))

# load syntelog table
syn <- as.data.frame(read_delim("../yucca_synteny/reciprocal_syntelogs_sameinYaYf_long.txt",delim = "\t"))

# make a function to format the parental data
format_parental_counts <- function(cdf){
  # where cdf = a parental counts data.frame
  ## remove metadata from cdf
  rownames(cdf) <- cdf$sample_name
  cdf <- cdf %>% select(-c(sample_name,genotype,time,treat,ZT,species))
  ## transpose cdf & make a GeneID column from the resultant row names
  tcdf <- as.data.frame(t(cdf))
  tcdf$GeneID <- rownames(tcdf)
  
  # merge with syntelog data to get a syntelogID; set this as the row names
  syncdf <- merge(tcdf,syn,by = "GeneID")
  rownames(syncdf) <- syncdf$syntelogID
  syncdf <- syncdf %>% select(-c(GeneID,syntelogID))
  return(syncdf)
}

ya_reform <- format_parental_counts(yacounts)
yf_reform <- format_parental_counts(yfcounts)

# merge data from all 3 species
## first: re-merge Yg samples (so there is a single counts value for each syntelog)
ygcountsonly <- ygcounts %>% select(-syntelogID)
part_strings <- unique(substr(colnames(ygcountsonly),1,nchar(colnames(ygcountsonly))-1))

ygcounts_sum <- data.frame(row.names = rownames(ygcountsonly))
for(i in part_strings){
  
  cna <- paste0(i, "a")
  cnf <- paste0(i, "f")
  
  has_a <- cna %in% colnames(ygcountsonly)
  has_f <- cnf %in% colnames(ygcountsonly)
  
  if(has_a & has_f){
    # Both exist → sum
    ygcounts_sum[[i]] <- ygcountsonly[[cna]] + ygcountsonly[[cnf]]
    
  } else if(has_a){
    # Only a exists → copy
    ygcounts_sum[[i]] <- ygcountsonly[[cna]]
    
  } else if(has_f){
    # Only f exists → copy
    ygcounts_sum[[i]] <- ygcountsonly[[cnf]]
  }
}

ygcounts_sum$syntelogID <- ygcounts$syntelogID
ya_reform$syntelogID <- rownames(ya_reform)
yf_reform$syntelogID <- rownames(yf_reform)
allcounts <- merge(ygcounts_sum,merge(ya_reform,yf_reform,by="syntelogID"),by="syntelogID")

# set syntelog IDs as row names & drop syntelogID column
rownames(allcounts) <- allcounts$syntelogID
allcounts <- allcounts %>% select(-syntelogID)

# set up the coldata object
## pull out metadata for each species
yacounts$genotype <- yacounts$species
yfcounts$genotype <- yfcounts$species
yamd <- yacounts %>% select(c(sample_name,genotype,time,treat,ZT,species))
yamd$generation <- "P1"
yamd$ploidy <- "diploid"
yfmd <- yfcounts %>% select(c(sample_name,genotype,time,treat,ZT,species))
yfmd$generation <- "P2"
yfmd$ploidy <- "diploid"

ygcmd <- as.data.frame(read_delim("./counts/Yg_toYgIS_all_correctedmd_countsmatrix.txt",delim = "\t"))
ygmd <- ygcmd %>% select(c(sample_name,genotype,time,treat,ZT,species))
ygmd$generation <- "F1"
ygmd$ploidy <- "homoploid"

# stick all the metadata together
allmd <- bind_rows(yamd,yfmd,ygmd)

# subset allmd to only the samples in allcounts
allmd <- allmd %>% filter(sample_name %in% colnames(allcounts))

# set sample_name as row names
rownames(allmd) <- allmd$sample_name
allmd <- allmd %>% select(-c(sample_name,time,species))

# create the SummarizedExperiment object
hybexp_se <- SummarizedExperiment(
  assays = list(counts = allcounts),
  colData = allmd[match(colnames(allcounts),rownames(allmd)),]
)

# save the SummarizedExperiment as an rds file