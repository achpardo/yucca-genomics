# Purpose: Run GO term enrichment for intersections of HybridExpress & time-structured genes.
# Author: Anna Pardo
# Date initiated: May 11, 2026

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

# combine to make Yg GO annotation
yggo <- c(yfgolist,yagolist)

# load lists of genes for enrichment
camhets <- fromJSON(file = "./CAM_TS_HybExp_classints.json")
c3hets <- fromJSON(file = "./C3+CAM_TS_HybExp_classints.json")
fachets <- fromJSON(file = "./facultative_CAM_TS_HybExp_classints.json")

format_genelist <- function(inlist,ptype){
  for(i in 1:length(names(inlist))){
    n <- names(inlist)[[i]]
    names(inlist)[[i]] <- paste0(ptype,"_TS_",n)
  }
  return(inlist)
}

camhets <- format_genelist(camhets,"CAM")
c3hets <- format_genelist(c3hets,"C3+CAM")
fachets <- format_genelist(fachets,"facultative_CAM")

hets <- c(camhets,c3hets,fachets)

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

hetsgo <- go_enrichment(hets,yggo)
unique(hetsgo$GO.term[!(hetsgo$GO.term %in% GO$id)])
hetsgo$Descr <- apply(hetsgo,MARGIN=1,function(x){
  get_term_property(GO,"name",x[2])})
write.table(hetsgo, file = "./GOenrichres_HybExp_TS_intgenes.txt",
            sep = "\t",row.names = FALSE,col.names = TRUE)
