# Purpose: Run GO enrichment on ASE & HybridExpress results (separately for Ya & Yf genes).
# Author: Anna Pardo
# Date initiated: May 1, 2026

# Load modules
library(topGO)
library(rjson)
library(dplyr)
library(ontologyIndex)
library(stringr)

setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/")
# load go.obo, which was downloaded May 1, 2026 from http://geneontology.org/docs/download-ontology/#go_obo_and_owl
GO <- get_ontology("./go.obo")

# load Ya & Yf GO annotations
yagolist <- fromJSON(file = "./Yav2_GO_annotation.json")
yfgolist <- fromJSON(file = "./Yfv3_hap1_GO_annotation.json")

# load ASE & HybridExpress results
