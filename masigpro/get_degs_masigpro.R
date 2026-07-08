# Purpose: Get drought vs. watered DEGs from maSigPro results (for each permutation).
# Author: Anna Pardo
# Date initiated: Apr. 21, 2026

# load modules
library(maSigPro)
library(readr)
library(dplyr)
library(edgeR)

setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/masigpro/C3+CAM_newsets_inputdata_maSigPro/")

# Apr. 17: load C3+CAM counts
ct11 <- as.data.frame(read_delim("./C3+CAM_11_3reps_msp_input.txt",delim = "\t"))
ct12 <- as.data.frame(read_delim("./C3+CAM_12_3reps_msp_input.txt",delim = "\t"))
ct13 <- as.data.frame(read_delim("./C3+CAM_13_3reps_msp_input.txt",delim = "\t"))

# make a function to format/normalize the data
make_tmm <- function(counts){
  # drop metadata columns
  countsonly <- counts %>% select(-any_of(c("genotype","time","treat","ZT","species","Total_Reads","phys")))
  countsonly <- as.data.frame(countsonly)
  row.names(countsonly) <- countsonly$sample_name
  countsonly <- countsonly %>% select(-sample_name)
  tcounts <- as.matrix(t(countsonly))
  dyg <- DGEList(tcounts)
  normdyg <- calcNormFactors(dyg,method = "TMM")
  normcounts <- cpm(normdyg,normalized.lib.sizes = TRUE)
  return(normcounts)
}

# function to construct watered & drought columns for metadata
treatmap <- function(t,md){
  maplist <- list()
  for(i in unique(md$treat)){
    if(i==t){
      maplist[i] <- 1
    } else {
      maplist[i] <- 0
    }
  }
  if(t=="W"){
    md$Watered <- as.integer(maplist[match(md$treat,names(maplist))])
  } else {
    md$Drought <- as.integer(maplist[match(md$treat,names(maplist))])
  }
  return(md)
}

# metadata setup function
metadata_setup <- function(counts){
  # create metadata object
  md <- counts %>% select(c("sample_name","genotype","treat","ZT"))
  # set sample name as row names
  md <- as.data.frame(md)
  row.names(md) <- md$sample_name
  # set ZT as Time column
  md$Time <- md$ZT
  # make a condition column
  md$condition <- paste0(md$treat,"_",md$ZT)
  # for each unique condition, assign a number & map this as a Replicate column
  repmap <- list()
  x = 1
  for(i in unique(md$condition)){
    repmap[i] <- x
    x = x+1
  }
  md$Replicate <- as.integer(repmap[match(md$condition,names(repmap))])
  md <- treatmap("W",md)
  md <- treatmap("D",md)
  
  # subset to relevant columns
  mspmd <- md %>% select(c(Time,Replicate,Watered,Drought))
  # make the design matrix
  design <- make.design.matrix(mspmd,degree = 5)
  return(design)
}

# function to run maSigPro
run_masigpro <- function(normcounts,design){
  fit <- p.vector(normcounts, design, Q = 0.05, MT.adjust = "BH",counts = TRUE,epsilon = 0.05)
  # find & remove influential genes
  tfit<-T.fit(data=fit,epsilon = 0.05)
  influential<-tfit$influ.info
  inf.genenames<-colnames(influential)
  nc2<-normcounts[!rownames(normcounts) %in% inf.genenames, ]
  # re-run p.vector()
  fit2 <- p.vector(nc2,design,Q=0.05,MT.adjust="BH",counts = TRUE,epsilon = 0.05)
  # and re-run T.fit()
  NBt <- T.fit(fit2,epsilon = 0.05)
  return(NBt)
}

run_all_de <- function(counts){
  tmm <- make_tmm(counts)
  des <- metadata_setup(counts)
  NBt <- run_masigpro(tmm,des)
  sigswd <- get.siggenes(NBt,rsq=0.7,vars="groups")
  return(sigswd)
}

setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/masigpro/DEGs_results_Jul06/")
degs13 <- run_all_de(ct13)
WDgenes<-degs13$sig.genes$DroughtvsWatered$sig.profiles
write.table(file="../DE_results/Jul06_degs_C3+CAM_13.txt", WDgenes)
WDgenesR2<-degs13$sig.genes$DroughtvsWatered$sig.pvalues$`R-squared`
write.table(file="../DE_results/Jul06_degs_C3+CAM_13_R2.txt", WDgenesR2)
WDgenespval<-degs13$sig.genes$DroughtvsWatered$sig.pvalues$`p-value`
write.table(file="../DE_results/Jul06_degs_C3+CAM_13_pval.txt", WDgenespval)
rm(degs13,WDgenes,WDgenesR2,WDgenespval)

degs12 <- run_all_de(ct12)
WDgenes<-degs12$sig.genes$DroughtvsWatered$sig.profiles
write.table(file="./Jul06_degs_C3+CAM_12.txt", WDgenes)
WDgenesR2<-degs12$sig.genes$DroughtvsWatered$sig.pvalues$`R-squared`
write.table(file="./Jul06_degs_C3+CAM_12_R2.txt", WDgenesR2)
WDgenespval<-degs12$sig.genes$DroughtvsWatered$sig.pvalues$`p-value`
write.table(file="./Jul06_degs_C3+CAM_12_pval.txt", WDgenespval)
rm(degs12,WDgenes,WDgenesR2,WDgenespval)

degs11 <- run_all_de(ct11)
WDgenes<-degs11$sig.genes$DroughtvsWatered$sig.profiles
write.table(file="../DE_results/Jul06_degs_C3+CAM_11.txt", WDgenes)
WDgenesR2<-degs11$sig.genes$DroughtvsWatered$sig.pvalues$`R-squared`
write.table(file="../DE_results/Jul06_degs_C3+CAM_11_R2.txt", WDgenesR2)
WDgenespval<-degs11$sig.genes$DroughtvsWatered$sig.pvalues$`p-value`
write.table(file="../DE_results/Jul06_degs_C3+CAM_11_pval.txt", WDgenespval)

# repeat for parental species
## load parental counts
setwd("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/counts/")
yacounts <- read_delim("./Ya_counts_withmd.txt",delim = "\t")
setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/masigpro/")
yfcounts <- read_delim("./inputdata_3reps/YfH1_counts_withmd_3reps_random.txt",delim = "\t")

setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/masigpro/C3+CAM_newsets_inputdata_maSigPro/")
degsyf <- run_all_de(yfcounts)
WDgenes<-degsyf$sig.genes$DroughtvsWatered$sig.profiles
write.table(file="../DE_results/degs_Yf.txt", WDgenes)
WDgenesR2<-degsyf$sig.genes$DroughtvsWatered$sig.pvalues$`R-squared`
write.table(file="../DE_results/degs_Yf_R2.txt", WDgenesR2)
WDgenespval<-degsyf$sig.genes$DroughtvsWatered$sig.pvalues$`p-value`
write.table(file="../DE_results/degs_Yf_pval.txt", WDgenespval)

degsya <- run_all_de(yacounts)
WDgenes<-degsya$sig.genes$DroughtvsWatered$sig.profiles
write.table(file="../DE_results/degs_Ya.txt", WDgenes)
WDgenesR2<-degsya$sig.genes$DroughtvsWatered$sig.pvalues$`R-squared`
write.table(file="../DE_results/degs_Ya_R2.txt", WDgenesR2)
WDgenespval<-degsya$sig.genes$DroughtvsWatered$sig.pvalues$`p-value`
write.table(file="../DE_results/degs_Ya_pval.txt", WDgenespval)
