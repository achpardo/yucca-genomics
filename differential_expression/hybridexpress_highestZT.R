# Purpose: Repeat HybridExpress analysis, but for each genotype/treatment combination (and parental), subset to a top
## percentage of highest expression.
# Author: Anna Pardo
# Date initiated: Mar. 12, 2026

library(HybridExpress)
library(DESeq2)
library(readr)
library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(ggplot2)
library(rjson)
library(matrixStats)
library(purrr)

setwd("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/")

# load summarizedExperiment object from previous HybridExpress run
se <- readRDS("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/differential_expression/summExp_forhybridExpress.rds")

## get vector of offspring genotypes
sp <- c("aloifolia","filamentosa","midparent")
hybgt <- unique(se$genotype)
hybgt <- hybgt[!hybgt %in% sp]

# function adjusted by ChatGPT
degs_treatgt <- function(gt, treatment, toppct = 20, expobj = se) {
  
  # Subset the SummarizedExperiment object
  gcd <- as.data.frame(colData(expobj)) %>%
    filter(genotype %in% c(gt, "aloifolia", "filamentosa")) %>%
    filter(treat == treatment)
  #print(unique(gcd$Generation))
  #print(gt)
  
  # Check if all generations are present
  if(!all(c("P1", "P2", "F1") %in% gcd$Generation)) {
    return("not all generation levels present")
  }
  
  # Prepare list to store top values for each genotype
  matl <- list()
  
  for(i in seq_along(unique(gcd$genotype))) {
    cur_gt <- unique(gcd$genotype)[i]
    
    subcd <- gcd %>% filter(genotype == cur_gt)
    samps <- rownames(subcd)
    gtassay <- assay(expobj)[, samps, drop = FALSE]
    
    # Number of top samples per row
    ntop <- round(ncol(gtassay) * (toppct / 100))
    if(ntop<3){
      ntop=3
    }
    print(paste0("Genotype/Species: ",cur_gt,", number of top samples: ",ntop))
    
    # Get top 'ntop' values per row
    topvals <- t(apply(gtassay, 1, function(x) sort(x, decreasing = TRUE)[1:ntop]))
    topvals <- as.data.frame(topvals)
    
    # Name columns: genotype_1, genotype_2, ...
    colnames(topvals) <- paste0(cur_gt, "_", seq_len(ntop))
    
    # Add syntelogID column for merging
    topvals$syntelogID <- rownames(topvals)
    
    matl[[i]] <- topvals
  }
  
  # Merge all genotype tables by syntelogID
  print("Merging all top values...")
  mergedAssay <- reduce(matl, inner_join, by = "syntelogID")
  
  # Remove syntelogID column from assay, keep for rownames
  rownames(mergedAssay) <- mergedAssay$syntelogID
  mergedAssay$syntelogID <- NULL
  
  # Set up new colData
  print("Setting up new colData...")
  prefix_map <- c("aloifolia" = "P1", "filamentosa" = "P2")
  prefix_map[gt] <- "F1"
  prefix <- sub("_.*", "", colnames(mergedAssay))
  print(prefix_map[prefix])
  
  # Safety check
  if(any(is.na(prefix_map[prefix]))) {
    stop("Some column names could not be mapped to Generation. Check column naming!")
  }
  
  cdata <- data.frame(Generation = prefix_map[prefix], genotype = prefix, row.names = colnames(mergedAssay))
  
  # Build SummarizedExperiment
  print("Building SummarizedExperiment...")
  gse <- SummarizedExperiment(assays = list(counts = as.matrix(mergedAssay)),
                              colData = cdata)
  
  # Add midparent expression and DESeq size factor normalization
  print("Adding midparent values & normalizing...")
  gse <- add_midparent_expression(gse)
  assay(gse) <- as.matrix(assay(gse))
  gse <- add_size_factors(gse)
  
  # Get differential expression list
  print("Getting DEGs...")
  deg_list <- get_deg_list(gse, alpha = 0.05)
  
  return(deg_list)
}

dgtdegs <- list()
wgtdegs <- list()
for(i in 1:length(hybgt)){
  dgtdegs[[i]] <- degs_treatgt(hybgt[i],"D",30)
  names(dgtdegs)[[i]] <- paste0(hybgt[i],"_D")
  
  wgtdegs[[i]] <- degs_treatgt(hybgt[i],"W",30)
  names(wgtdegs)[[i]] <- paste0(hybgt[i],"_W")
}

# check which genotype_treatments didn't work for drought
remove <- c()
for(i in 1:length(dgtdegs)){
  if(is.character(dgtdegs[[i]])){
    print(names(dgtdegs)[[i]])
    remove <- append(remove,names(dgtdegs)[[i]])
  }
  if(is.character(wgtdegs[[i]])){
    print(names(wgtdegs)[[i]])
    remove <- append(remove,names(wgtdegs)[[i]])
  }
}

# make a better list object
alldegs <- append(dgtdegs,wgtdegs)
alldegs <- alldegs[names(alldegs) %in% remove == FALSE]

# create summary information
allsums <- list()
for(i in 1:length(alldegs)){
  allsums[[i]] <- get_deg_counts(alldegs[[i]])
  allsums[[i]]$genotype_treat <- names(alldegs)[[i]]
  names(allsums)[[i]] <- names(alldegs)[[i]]
}
allsum <- bind_rows(allsums)
write_delim(allsum,file = "//wsl$/Ubuntu/home/leviathan22/yucca-genomics/differential_expression/summary_hybexp_top30pct_degs_bygttreat.txt",delim = "\t")

# do expression partitioning for each DEG set
allparts <- list()
for(i in 1:length(alldegs)){
  allparts[[i]] <- expression_partitioning(alldegs[[i]])
  allparts[[i]]$genotype_treat <- names(alldegs)[[i]]
}
apdf <- bind_rows(allparts)
