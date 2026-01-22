# Purpose: Convert all counts generated from unique alignments to Y. gloriosa in silico genome to TPM.
# Author: Anna Pardo
# Date initiated: Sept. 15, 2025

# set working directory
setwd("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/YgIS_counts/")

# load modules
library(edgeR)
library(dplyr)
library(readr)
library(DGEobj.utils)
library(janitor)
library(purrr)

# load data
# load counts matrices: Yg mapped to Yg in silico genome
files <- list.files(
  path = "//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/YgIS_counts/",
  pattern = "\\.txt$",
  full.names = TRUE
)
dflist <- list()

b = 1
for(f in files){
  data <- read_delim(f,delim = "\t",comment = "#")
  dflist[[b]] <- data %>% select(-c(Chr,Start,End,Strand))
  b = b+1
}

# merge all files together
## start with the first two
ygmerged <- merge(dflist[[1]],dflist[[2]])
for(i in c(3:38)){
  ygmerged <- merge(ygmerged,dflist[[i]])
}

# rename the columns
newcolnames <- vector()
for (i in colnames(ygmerged)) {
  if(i=="Geneid"){
    newcolnames <- append(newcolnames,"GeneID")
  } else if(i=="Length") {
    newcolnames <- append(newcolnames,"Length")
  } else {
    x <- basename(i)
    y <- strsplit(x,"_u")[[1]][1]
    z <- strsplit(y,"_p")[[1]][1]
    if(startsWith(z,"SRR")){
      z <- strsplit(z,"_")[[1]][1]
    }
    newcolnames <- append(newcolnames,z)
  }
}
colnames(ygmerged) <- newcolnames

# now split off the counts matrix for convertCounts() (all columns must be numeric)
ygcm <- ygmerged %>% select(-c(GeneID,Length))

# convert all columns to numeric
ygcmNum <- ygcm %>% mutate_all(function(x) as.numeric(as.character(x)))

# check for NAs
ygcmNum[which(is.na(ygcmNum)),]
# none present

# check for columns with a sum of 0
cs <- colSums(ygcmNum)
x <- which(cs==0)
x
# Y318 has no counts - drop this sample
ygcmNum <- ygcmNum %>% select(!(Y318))

# and split off a new geneLengths object
geneLengths_yg <- ygmerged %>% select(c(GeneID,Length))

# convert counts to TPM
ygtpm <- convertCounts(countsMatrix = as.matrix(ygcmNum),unit = "tpm",geneLength = geneLengths_yg$Length,log = FALSE,normalize = "none")

# add back gene IDs and save as txt file 
ygtpm <- data.frame(ygtpm)
ygtpm$GeneID <- ygmerged$GeneID
write_delim(ygtpm,file = "//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/TPM/Yg_toYgIS_allTPM_nomd_15-Sep-2025.txt",delim = "\t")
# save also: ygmerged
write_delim(ygmerged,file = "//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/counts/Yg_toYgIS_all_nomd_countsmatrix_15-Sep-2025.txt",delim = "\t")
