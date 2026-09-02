# Purpose: Run GO term enrichment on TSGs shared & unique among Y. gloriosa physiotypes.
# Author: Anna Pardo
# Date initiated: Sept. 1, 2026

# Load modules
library(topGO)
library(rjson)
library(dplyr)
library(ontologyIndex)
library(stringr)
library(readr)

setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/")
# load go.obo, which was downloaded May 1, 2026 from http://geneontology.org/docs/download-ontology/#go_obo_and_owl
GO <- get_ontology("./go.obo")

# load parental GO annotations
yago <- fromJSON(file = "./Yav2_GO_annotation.json")
yfgo <- fromJSON(file = "./Yfv3_hap1_GO_annotation.json")
# make combined GO annotation - Yg overall
yggo <- c(yfgo,yago)

# load gene lists: 3-way shared & unique to each physiotype TSGs
tsgs <- read_delim("./masigpro/3way_TSGs_Yg_allphys.txt",delim = "\t",col_names = F)$X1
ucam <- read_delim("./masigpro/unique_CAM_TSGs.txt",delim = "\t",col_names = F)$X1
uc3 <- read_delim("./masigpro/unique_C3+CAM_TSGs.txt",delim = "\t",col_names = F)$X1
ufac <- read_delim("./masigpro/unique_facCAM_TSGs.txt",delim = "\t",col_names = F)$X1

# need to run GO enrichment on these gene sets, excluding the pathway genes. to do so, load pathway annotation
pathgenes <- as.data.frame(read_csv(file = "./photresp_N_citrate_genes_Yucca.csv"))
clock <- as.data.frame(read_csv("./circadian_light_genes_by_orthology_Yucca.csv"))
cam <- as.data.frame(read_delim("./camgenes_Ya_Yf_orthology_synteny.txt",delim = "\t"))
# put together a list of all the genes in these pathways (vector)
pgenes_vec <- unique(c(clock$GeneID,cam$GeneID,pathgenes$GeneID))

# subset TSG sets
tsgs_nopath <- setdiff(tsgs,pgenes_vec)
ucam_nopath <- setdiff(ucam,pgenes_vec)
ufac_nopath <- setdiff(ufac,pgenes_vec)
uc3_nopath <- setdiff(uc3,pgenes_vec)

# combine these into a named list of gene sets
genesets <- list(tsgs_nopath,ucam_nopath,ufac_nopath,uc3_nopath)
names(genesets) <- c("3-way_overlap","CAM_unique","facCAM_unique","C3+CAM_unique")

# GO enrichment function #####
go_enrichment <- function(injson,spgo){
  ### injson is input json file (named list) with genes of interest ###
  ### spgo is the GO annotation for the species (as named list) ###
  
  allgenes_sp <- names(spgo)
  intgenes <- lapply(injson,function(x){factor(as.integer(allgenes_sp %in% x))})
  
  valid <- sapply(intgenes, function(x) length(unique(x)) == 2)
  intgenes <- intgenes[valid]
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
  all.sig.GO.bp <- all.sig.GO.bp[lapply(all.sig.GO.bp,nrow)>0]
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

# run GO enrichment
all_enrich <- go_enrichment(genesets,yggo)
# now add the descriptions of the GO terms
all_enrich$Descr <- apply(all_enrich,MARGIN=1,function(x){
  get_term_property(GO,"name",x[2])})
# save final results as .txt file
write.table(all_enrich, file = "./masigpro/masigpro_TS_shared_unique_GOenrich.txt",
            sep = "\t",row.names = FALSE,col.names = TRUE)
