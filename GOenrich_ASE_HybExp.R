# Purpose: Run GO enrichment on ASE & HybridExpress results (separately for Ya & Yf genes).
# Author: Anna Pardo
# Date initiated: May 1, 2026

# Load modules
library(topGO)
library(rjson)
library(dplyr)
library(ontologyIndex)
library(stringr)
library(tidyr)

setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/")
# load go.obo, which was downloaded May 1, 2026 from http://geneontology.org/docs/download-ontology/#go_obo_and_owl
GO <- get_ontology("./go.obo")

# load Ya & Yf GO annotations
yagolist <- fromJSON(file = "./Yav2_GO_annotation.json")
yfgolist <- fromJSON(file = "./Yfv3_hap1_GO_annotation.json")

# load ASE & HybridExpress results
ase <- read.csv("./differential_expression/ASE_gttreat_output.txt",sep = "\t")
hybexp <- read.csv("./differential_expression/expression_partitioning_downstream/hybexp_top30pct_partitiongenes_bygttreat.txt",
                   sep = "\t")

# for each of these: split out genotype & treatment information
hybexp <- hybexp %>% separate_wider_delim(genotype_treat,delim = "_",names = c("genotype","treat"))
ase <- ase %>% separate_wider_delim(Contrast,delim = "-V-",names = c("ref","cont"))
ase <- ase %>% separate_wider_delim(ref,delim = "_",names = c("genotype","treat","sp")) %>% select(-c("sp","cont"))

# load physiology information
phys <- fromJSON(file = "./physiology/physiological_categories_from_TA.json")

# add physiology columns to hybexp & ase
ase <- ase %>% mutate(CAMphys = recode(genotype,!!!phys))
hybexp <- hybexp %>% mutate(CAMphys = recode(genotype,!!!phys))

# for ASE: if log2 fold change is negative, gene is biased to Ya; if positive, biased to Yf
ase <- ase %>% mutate(bias = if_else(log2FoldChange > 0, "Yf", "Ya"))

# combine phys & treatment info
ase$phys_treat <- paste0(ase$CAMphys,"_",ase$treat)
hybexp$phys_treat <- paste0(hybexp$CAMphys,"_",hybexp$treat)

# load syntelog information
synwide <- read.delim("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/yucca_synteny/reciprocal_syntelogs_sameinYaYf.txt",sep="\t")
synlong <- read.delim("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/yucca_synteny/reciprocal_syntelogs_sameinYaYf_long.txt",sep = "\t")

# set up the gene sets list: by physiology-treatment set
hybexp$class_phys_treat <- paste0(hybexp$phys_treat,"_",hybexp$Class)
cpt <- unique(hybexp$class_phys_treat)
he_genes_list <- list()
for(i in 1:length(cpt)){
  df <- hybexp %>% filter(class_phys_treat==cpt[i])
  he_genes_list[[i]] <- unique(df$Gene)
  names(he_genes_list)[[i]] <- cpt[i]
} 

# for HybridExpress: convert GO annotations to syntelog-wise, i.e. include all GO annotations for both syntelogs
yasyn <- list()
yfsyn <- list()
for (i in 1:length(synwide$syntelogID)) {
  yasyn[[i]] <- synwide$Ya_GeneID[i]
  yfsyn[[i]] <- synwide$Yf_GeneID[i]
  names(yasyn)[[i]] <- synwide$syntelogID[i]
  names(yfsyn)[[i]] <- synwide$syntelogID[i]
}


# collapse duplicates by gene
yfgolist <- tapply(unlist(yfgolist),
                   rep(names(yfgolist), lengths(yfgolist)),
                   unique,
                   simplify = FALSE)

common_ids <- intersect(names(yasyn), names(yfsyn))

syntelog_go <- setNames(
  lapply(common_ids, function(id) {
    unique(c(
      yagolist[[ yasyn[[id]] ]],
      yfgolist[[ yfsyn[[id]] ]]
    ))
  }),
  common_ids
)
syntelog_go <- syntelog_go[sapply(syntelog_go, length) > 0]

# save as JSON file
write(toJSON(syntelog_go),"./GOannotation_by_syntelog_ID.json")

# for ASE: set up Ya & Yf-biased gene sets by phys_treat
ase_genesets <- list()
ase$bias_phys_treat <- paste0(ase$phys_treat,"_",ase$bias)
ase$syntelogID <- ase$GeneID
ase <- merge(ase,synwide)
asegenes <- list()
bpt <- unique(ase$bias_phys_treat)
for(i in 1:length(bpt)){
  df <- ase %>% filter(bias_phys_treat==bpt[i])
  if(grepl("Yf",bpt[i])){
    asegenes[[i]] <- unique(df$Yf_GeneID)
  } else {
    asegenes[[i]] <- unique(df$Ya_GeneID)
  }
  names(asegenes)[[i]] <- bpt[i]
}

# GO enrichment function #####
go_enrichment <- function(injson,spgo){
  ### injson is input json file (named list) with genes of interest ###
  ### spgo is the GO annotation for the species (as named list) ###
  
  allgenes_sp <- names(spgo)
  intgenes <- lapply(injson,function(x){factor(as.integer(allgenes_sp %in% x))})
  
  valid <- sapply(intgenes, function(x) length(unique(x)) == 2)
  intgenes <- intgenes[valid]
  names(intgenes) <- names(injson)[valid]
  print(names(intgenes))
  
  for (i in 1:length(intgenes)){
    names(intgenes[[i]]) <- allgenes_sp
  }
  #print(str(intgenes))
  
  all.GOdata.bp <- lapply(intgenes,function(x){
    new("topGOdata",ontology="BP",allGenes=x,annot=annFUN.gene2GO,gene2GO=spgo)})
  
  # initialize test stat - Fisher's exact test
  test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test")
  
  all.resultFisher.bp <- lapply(all.GOdata.bp,function(x){getSigGroups(x,test.stat)})
  all.sig.GO.bp <- lapply(all.resultFisher.bp,function(x){data.frame("p-value"=score(x))})
  for (i in 1:length(all.sig.GO.bp)){
    all.sig.GO.bp[[i]]$GO.term <- row.names(all.sig.GO.bp[[i]])
    all.sig.GO.bp[[i]]$p.adj <- p.adjust(all.sig.GO.bp[[i]]$p.value,method = "fdr")
    all.sig.GO.bp[[i]] <- all.sig.GO.bp[[i]][all.sig.GO.bp[[i]]$p.adj<0.05,]
  }
  all.sig.GO.bp <- all.sig.GO.bp[sapply(all.sig.GO.bp,nrow)>0]
  for (i in 1:length(all.sig.GO.bp)){
    all.sig.GO.bp[[i]]$regset <- names(all.sig.GO.bp)[i]
    all.sig.GO.bp[[i]]$Ontology <- "BP"
  }
  for(i in 1:length(all.sig.GO.bp)){
    all.sig.GO.bp[[i]] %>% filter(p.adj < 0.05)
  }
  prelimdf <- do.call(rbind,all.sig.GO.bp)
  
  all.sig.GO.res.bp <- lapply(seq_along(all.resultFisher.bp),function(x,y,z){
    GenTable(z[[x]],y[[x]],topNodes = 250)},y=all.resultFisher.bp,z=all.GOdata.bp)
  #add module classification to each df in this list
  for(i in seq_along(all.sig.GO.res.bp)){
    all.sig.GO.res.bp[[i]]$regset <- names(all.resultFisher.bp[i])
  }
  # convert list to df
  sigann_df <- bind_rows(all.sig.GO.res.bp)
  # join to previous df
  dfname <- left_join(prelimdf,sigann_df,by = c("GO.term" = "GO.ID","regset"))
  
  return(dfname)
}

# run HybridExpress GO enrichment
hego <- go_enrichment(he_genes_list,syntelog_go)
unique(hego$GO.term[!(hego$GO.term %in% GO$id)])
# now add the descriptions of the GO terms
hego$Descr <- apply(hego,MARGIN=1,function(x){
  get_term_property(GO,"name",x[2])})
# save final results as .txt file
write.table(hego, file = "./differential_expression/expression_partitioning_downstream/HybExp_GOenrich_res.txt",
            sep = "\t",row.names = FALSE,col.names = TRUE)

# run ASE GO enrichment
yggo <- c(yfgolist,yagolist)
asego <- go_enrichment(asegenes,yggo)
unique(asego$GO.term[!(asego$GO.term %in% GO$id)])
asego$Descr <- apply(asego,MARGIN=1,function(x){
  get_term_property(GO,"name",x[2])})
write.table(asego, file = "./differential_expression/ASE_GOenrich_res.txt",
            sep = "\t",row.names = FALSE,col.names = TRUE)
