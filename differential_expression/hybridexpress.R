# Purpose: Try out HybridExpress for gene categorization.
# Author: Anna Pardo
# Date initiated: Feb. 11, 2026

library(HybridExpress)
library(DESeq2)
library(readr)
library(dplyr)
library(SummarizedExperiment)
library(ggplot2)
library(rjson)

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
yamd$Generation <- "P1"
yamd$Ploidy <- "diploid"
yfmd <- yfcounts %>% select(c(sample_name,genotype,time,treat,ZT,species))
yfmd$Generation <- "P2"
yfmd$Ploidy <- "diploid"

ygcmd <- as.data.frame(read_delim("./counts/Yg_toYgIS_all_correctedmd_countsmatrix.txt",delim = "\t"))
ygmd <- ygcmd %>% select(c(sample_name,genotype,time,treat,ZT,species))
ygmd$Generation <- "F1"
ygmd$Ploidy <- "homoploid"

# stick all the metadata together
allmd <- bind_rows(yamd,yfmd,ygmd)

# subset allmd to only the samples in allcounts
allmd <- allmd %>% filter(sample_name %in% colnames(allcounts))
allmd <- allmd %>% filter(sample_name %in% colnames(assay(hybexp_se)))

# set sample_name as row names
rownames(allmd) <- allmd$sample_name
allmd <- allmd %>% select(-c(sample_name,time,species))

# create the SummarizedExperiment object
hybexp_se <- SummarizedExperiment(
  assays = list(counts = allcounts),
  colData = allmd[match(colnames(allcounts),rownames(allmd)),]
)

hybexp_se2 <- SummarizedExperiment(
  assays = as.matrix(assay(hybexp_se)[1]),
  colData = allmd[match(colnames(assay(hybexp_se)),rownames(allmd)),]
)

# save the SummarizedExperiment as an rds file
saveRDS(hybexp_se2,"//wsl$/Ubuntu/home/leviathan22/yucca-genomics/differential_expression/summExp_forhybridExpress.rds")

######## Run HybridExpress: get DEGs #####
# add midparent values
hese <- add_midparent_expression(hybexp_se2)
assay(hese) <- as.matrix(assay(hese))
hese <- add_size_factors(hese)

# do some exploratory plotting with HybridExpress
# For colData rows with missing values (midparent samples), add "midparent"
hese$Ploidy[is.na(hese$Ploidy)] <- "midparent"
hese$Generation[is.na(hese$Generation)] <- "midparent"
hese$genotype[is.na(hese$genotype)] <- "midparent"
hese$treat[is.na(hese$treat)] <- "midparent"
hese$ZT[is.na(hese$ZT)] <- "midparent"

# custom color palette from https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

pca_plot(hese, color_by = "Generation", shape_by = "Ploidy", add_mean = TRUE)
pca_plot(hese, ntop=10000, color_by = "genotype", shape_by = "Generation", add_mean = TRUE,
         palette = c25) +
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10))

pca_plot(hese, ntop=10000, color_by = "genotype", shape_by = "treat", add_mean = TRUE,
         palette = c25) +
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10))

pca_plot(hese, ntop=10000, color_by = "ZT", shape_by = "Generation", add_mean = TRUE,
         palette = c25) +
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10))

# get DEGs: first (most basic) by generation, without additional metadata
de_bygen <- get_deg_list(hese,alpha = 0.05)

# run by genotype only
## get vector of offspring genotypes
sp <- c("aloifolia","filamentosa","midparent")
hybgt <- unique(hese$genotype)
hybgt <- hybgt[!hybgt %in% sp]
de_bygen_gt <- get_deg_list(hese,coldata_column = "genotype",parent1 = "aloifolia",
                            parent2 = "filamentosa",offspring = hybgt,
                            midparent = "midparent",alpha = 0.05)

# I don't think that did what I wanted it to...let's save both of these objects anyway:
saveRDS(de_bygen,"//wsl$/Ubuntu/home/leviathan22/yucca-genomics/differential_expression/hybexp_degenes_bygenerationonly.rds")
saveRDS(de_bygen,"//wsl$/Ubuntu/home/leviathan22/yucca-genomics/differential_expression/hybexp_degenes_bygeneration_gtcol_allgts.rds")

# set up a function to run this for each genotype individually
degs_genotype <- function(gt){
  # subset the summarizedexperiment object
  gcd <- as.data.frame(hese@colData) %>% filter(genotype %in% c(gt,"aloifolia","filamentosa","midparent"))
  print(head(gcd))
  samps <- rownames(gcd)
  gtassay <- assay(hese)[,samps]
  
  gse <- SummarizedExperiment(assays = gtassay,colData = gcd)
  
  deg_list <- get_deg_list(gse,alpha = 0.05)
  return(deg_list)
}

# loop through genotypes & get all lists
gtdegs <- list()
for(i in 1:length(hybgt)){
  gtdegs[[i]] <- degs_genotype(hybgt[i])
  names(gtdegs)[[i]] <- hybgt[i]
}

# save to an RDS just in case of session loss
saveRDS(gtdegs,"//wsl$/Ubuntu/home/leviathan22/yucca-genomics/differential_expression/hybexp_degenes_eachgt_bygen.rds")

# get DEG summaries for each genotype
gtsums <- list()
for(i in 1:length(gtdegs)){
  gtsums[[i]] <- get_deg_counts(gtdegs[[i]])
  names(gtsums)[[i]] <- names(gtdegs)[[i]]
}

# save as JSON