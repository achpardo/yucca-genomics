# Purpose: Run maSigPro with selected genotypes of Yg combined (either 1 or 2 reps per condition).
# Author: Anna Pardo
# Date initiated: Apr. 13, 2025

# load modules
library(maSigPro)
library(readr)
library(dplyr)
library(edgeR)

setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/masigpro/C3+CAM_newsets_inputdata_maSigPro/")
# load counts dataframes
c3ct <- as.data.frame(read_delim("./Yg_C3_counts_withmd_4reps_random.txt",delim="\t"))
pfcct <- as.data.frame(read_delim("./Yg_possible_fac_CAM_counts_withmd_4reps_random.txt",delim="\t"))
fcct <- as.data.frame(read_delim("./Yg_facultative_CAM_counts_withmd_4reps_random.txt",delim="\t"))

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

# function for setting up metadata
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
  
  # make a function to construct the Watered & Drought cols in the same way
  treatmap <- function(t){
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
  
  md <- treatmap("W")
  md <- treatmap("D")
  
  # subset to relevant columns
  mspmd <- md %>% select(c(Time,Replicate,Watered,Drought))
  # make the design matrix
  design <- make.design.matrix(mspmd,degree = 5)
  # degree = (number of unique time points) - 1
  return(design)
}

# function to run maSigPro
run_masigpro <- function(normcounts,design){
  fit <- p.vector(normcounts, design, Q = 0.05, MT.adjust = "BH")
  # find & remove influential genes
  tfit<-T.fit(data=fit)
  influential<-tfit$influ.info
  inf.genenames<-colnames(influential)
  nc2<-normcounts[!rownames(normcounts) %in% inf.genenames, ]
  # re-run p.vector()
  fit2 <- p.vector(nc2,design,Q=0.05,MT.adjust="BH")
  return(fit2)
}

# make wrapper function
run_all <- function(counts){
  tmm <- make_tmm(counts)
  des <- metadata_setup(counts)
  fit2 <- run_masigpro(tmm,des)
  prof <- see.genes.kh(fit2$SELEC, edesign=fit2$edesign, dis=fit2$dis, cluster.data=1, 
                       groups.vector=fit2$groups.vector, cluster.method="hclust", k=6)
  # extract genes in clusters
  clust <- as.data.frame(prof$cut)
  clust$GeneID <- row.names(clust)
  clust$cluster <- clust$`prof$cut`
  clust <- clust %>% select(c(cluster,GeneID))
  return(clust)
}

pfcclust <- run_all(pfcct)
fcclust <- run_all(fcct)
c3clust <- run_all(c3ct)

write.table(pfcclust,file = "../Yg_possible_fac_CAM_4reps_masigpro_clusters.txt",sep = "\t",col.names=T,row.names=F)
write.table(fcclust,file = "../Yg_facultative_CAM_4reps_masigpro_clusters.txt",sep = "\t",col.names=T,row.names=F)
write.table(c3clust,file = "../Yg_C3_4reps_masigpro_clusters.txt",sep = "\t",col.names=T,row.names=F)

# repeat for counts with 3 replicates
setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/masigpro/inputdata_3reps/")
yfct <- as.data.frame(read_delim("./YfH1_counts_withmd_3reps_random.txt",delim = "\t"))
c3ct <- as.data.frame(read_delim("./Yg_C3_counts_withmd_3reps_random.txt",delim = "\t"))
pfcct <- as.data.frame(read_delim("./Yg_possible_fac_CAM_counts_withmd_3reps_random.txt",delim = "\t"))
fcct <- as.data.frame(read_delim("./Yg_facultative_CAM_counts_withmd_3reps_random.txt",delim = "\t"))

yfclust <- run_all(yfct)
write.table(yfclust,file = "../Yf_3reps_masigpro_clusters.txt",sep = "\t",col.names = T,row.names = F)

c3clust <- run_all(c3ct)
write.table(c3clust,file = "../Yg_C3_3reps_masigpro_clusters.txt",sep = "\t",col.names = T,row.names = F)

pfcclust <- run_all(pfcct)
write.table(pfcclust,file = "../Yg_possible_fac_CAM_3reps_masigpro_clusters.txt",sep = "\t",col.names = T,row.names = F)

fcclust <- run_all(fcct)
write.table(fcclust,file = "../Yg_facultative_CAM_3reps_masigpro_clusters.txt",sep = "\t",col.names = T,row.names = F)

# April 17: Run for C3+CAM counts subsets
clust11 <- run_all(ct11)
write.table(clust11,file = "./C3+CAM_11_3reps_msp_input.txt_clusters.txt",sep = "\t",col.names = T,row.names = F)

clust12 <- run_all(ct12)
write.table(clust12,file = "./C3+CAM_12_3reps_msp_input.txt_clusters.txt",sep = "\t",col.names = T,row.names = F)

clust13 <- run_all(ct13)
write.table(clust13,file = "./C3+CAM_13_3reps_msp_input.txt_clusters.txt",sep = "\t",col.names = T,row.names = F)
