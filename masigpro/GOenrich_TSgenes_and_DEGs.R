# Purpose: Run GO term enrichment for time-structured & DE gene sets from maSigPro.
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

# load annotations for Ya & Yf genomes
yfannot <- read.delim("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/gene.functions.txt",sep = "\t")
# keep only IdType=="GO"
yfgo <- yfannot %>% filter(IdType=="GO")
yfgo$Description <- yfgo$description.if.any..GO.is.via.pfam2go.interpro
yfgo <- yfgo %>% select(c(GeneID,Id,IdType,Description))
# write this to a file
write.table(yfgo,file = "./Yfv3_hap1_GO_annotation.txt",sep = "\t",col.names = T,row.names = F)

# load file containing GO annotation for Ya (and subset to relevant columns)
yaannot <- read.csv("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/Y_aloifolia.emapper.annotations.csv",sep = ",")
yago <- yaannot %>% select(c(query,GOs))
#yago$query <- gsub("\\s+", "",yago$query)
# gene IDs are in the wrong format (query column) - convert to proper ending
yago$GeneID <- gsub("\\.1\\.p$", ".v2.1", yago$query)
## filter yago to only those genes which actually have GO annotations
yago <- yago %>% filter(GOs!="-")
# set up: named list where name=GeneID, value=character vector of GO term IDs
yagolist <- list()
for(i in 1:length(row.names(yago))){
  s <- yago$GOs[[i]] # this is a single string containing all GO IDs
  yagolist[[i]] <- str_split_1(s,",")
  names(yagolist)[[i]] <- yago$GeneID[i]
}
# write to JSON
write(toJSON(yagolist),"./Yav2_GO_annotation.json")

yfgolist <- list()
for(i in 1:length(unique(yfgo$GeneID))){
  df <- yfgo %>% filter(GeneID==yfgo$GeneID[i])
  yfgolist[[i]] <- as.character(df$Id)
  names(yfgolist)[[i]] <- yfgo$GeneID[i]
}
write(toJSON(yfgolist),"./Yfv3_hap1_GO_annotation.json")

# make combined GO annotation - Yg overall
yggo <- c(yfgolist,yagolist)

# load data
## time-structured gene sets
# final gene set format needed: named list where name=gene set name, value=character vector of gene IDs
camts <- read.csv("./masigpro/CAM_TSgenes_9gt.txt",sep = "\t")
c3ts <- read.csv("./masigpro/C3+CAM_TSgenes_9gt.txt",sep = "\t")
facts <- read.csv("./masigpro/facCAM_TSgenes_9gt.txt",sep = "\t")
yfts <- read.csv("./masigpro/Yf_3reps_masigpro_clusters.txt",sep = "\t")
yats <- read.csv("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/masigpro_results/Ya_hclust_k6_clusters.txt",
                 sep = "\t")

# load DEG sets
yade <- as.character(read.csv("./masigpro/Ya_DEGs_list.txt",sep = "\t",header = F)$V1)
yfde <- as.character(read.csv("./masigpro/Yf_DEGs_list.txt",sep = "\t",header = F)$V1)
ygde <- read.csv("./masigpro/Yg_maSigPro_DEGs_allphys.txt",sep = "\t")
camde <- ygde %>% filter(phys=="CAM")
c3de <- ygde %>% filter(phys=="C3+CAM")
facde <- ygde %>% filter(phys=="facultative_CAM")

genesets <- list(camts$GeneID,c3ts$GeneID,facts$GeneID,yfts$GeneID,yats$GeneID,
                 yade,yfde,camde$GeneID,c3de$GeneID,facde$GeneID)
names(genesets) <- c(
  "YgCAM_TS",
  "YgC3+CAM_TS",
  "YgfacCAM_TS",
  "Yf_TS",
  "Ya_TS",
  "Ya_DE",
  "Yf_DE",
  "YgCAM_DE",
  "YgC3+CAM_DE",
  "YgfacCAM_DE"
)

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
# check to see if there are any non-matching GO terms 
unique(all_enrich$GO.term[!(all_enrich$GO.term %in% GO$id)])
# no mismatches present

# now add the descriptions of the GO terms
all_enrich$Descr <- apply(all_enrich,MARGIN=1,function(x){
  get_term_property(GO,"name",x[2])})

# save final results as .txt file
write.table(all_enrich, file = "./masigpro/masigpro_TS_DE_bothsg_goenrich_01-May-2026.txt",
            sep = "\t",row.names = FALSE,col.names = TRUE)
