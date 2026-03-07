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
