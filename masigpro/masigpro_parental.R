# Purpose: Run maSigPro by treatment on both parental species (Ya and Yf).
# Author: Anna Pardo
# Date initiated: Nov. 10, 2025

# load modules
library(maSigPro)
library(readr)
library(dplyr)
library(edgeR)
library(Mfuzz)

######### Y. aloifolia ##########
setwd("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/counts/")
mdcounts <- read_delim("./Ya_counts_withmd.txt",delim = "\t")

# normalize counts with edgeR
# first: format counts for this normalization
# check for columns that are not gene IDs
for(c in colnames(mdcounts)){
  if(startsWith(c,"Yu")==FALSE){
    print(c)
  }
}
# drop metadata columns
counts <- mdcounts %>% select(-c("genotype","time","treat","ZT","species"))
counts <- as.data.frame(counts)
row.names(counts) <- counts$sample_name
counts <- counts %>% select(-sample_name)
tcounts <- as.matrix(t(counts))

dyg <- DGEList(tcounts)
normdyg <- calcNormFactors(dyg,method = "TMM")
normcounts <- cpm(normdyg,normalized.lib.sizes = TRUE)

# create metadata object
md <- mdcounts %>% select(c("sample_name","genotype","treat","ZT"))
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

fit <- p.vector(normcounts, design, Q = 0.05, MT.adjust = "BH")

# find & remove influential genes
tfit<-T.fit(data=fit)
influential<-tfit$influ.info
inf.genenames<-colnames(influential)
nc2<-normcounts[!rownames(normcounts) %in% inf.genenames, ]

# re-run p.vector()
fit2 <- p.vector(nc2,design,Q=0.05,MT.adjust="BH")

# pick k for k-means clustering (I think it will be 7 but let's verify...)
wss<-(nrow(fit2$SELEC)-1)*sum(apply(fit2$SELEC,2,var))
for (i in 2:15){
  wss[i]<- sum(kmeans(fit2$SELEC, centers=i, iter.max=20)$withinss)
} 
plot(1:15, wss, type="b")

pk7kmeans<-see.genes.kh(fit2$SELEC, edesign=fit2$edesign, dis=fit2$dis, cluster.data=1, groups.vector=fit2$groups.vector, cluster.method="kmeans", k=7)

phck7 <- see.genes.kh(fit2$SELEC, edesign=fit2$edesign, dis=fit2$dis, cluster.data=1, groups.vector=fit2$groups.vector, cluster.method="hclust", k=7)
phck5 <- see.genes.kh(fit2$SELEC, edesign=fit2$edesign, dis=fit2$dis, cluster.data=1, groups.vector=fit2$groups.vector, cluster.method="hclust", k=5)
pk5kmeans<-see.genes.kh(fit2$SELEC, edesign=fit2$edesign, dis=fit2$dis, cluster.data=1, groups.vector=fit2$groups.vector, cluster.method="kmeans", k=5)

phck6 <- see.genes.kh(fit2$SELEC, edesign=fit2$edesign, dis=fit2$dis, cluster.data=1, groups.vector=fit2$groups.vector, cluster.method="hclust", k=6)

# function for extracting & saving cluster data
saveclust <- function(profilesall,outfile){
  clust <- as.data.frame(profilesall$cut)
  clust$GeneID <- row.names(clust)
  clust <- clust %>% rename(cluster = `profilesall$cut`)
  write.table(clust,file = outfile,sep = "\t",col.names = T,row.names = F)
}

# save all 4 cluster objects
saveclust(pk7kmeans,"../masigpro_results/Ya_kmeans_k7_clusters.txt")
saveclust(pk5kmeans,"../masigpro_results/Ya_kmeans_k5_clusters.txt")
saveclust(phck5,"../masigpro_results/Ya_hclust_k5_clusters.txt")
saveclust(phck7,"../masigpro_results/Ya_hclust_k7_clusters.txt")

saveclust(phck6,"../masigpro_results/Ya_hclust_k6_clusters.txt")

############## Y. filamentosa #########
# load data
mdcounts <- read_delim("./YfH1_counts_withmd.txt",delim = "\t")

# normalize counts with edgeR
# first: format counts for this normalization
# check for columns that are not gene IDs
for(c in colnames(mdcounts)){
  if(startsWith(c,"Yu")==FALSE){
    print(c)
  }
}
# drop metadata columns
counts <- mdcounts %>% select(-c("genotype","time","treat","ZT","species"))
counts <- as.data.frame(counts)
row.names(counts) <- counts$sample_name
counts <- counts %>% select(-sample_name)
tcounts <- as.matrix(t(counts))

dyg <- DGEList(tcounts)
normdyg <- calcNormFactors(dyg,method = "TMM")
normcounts <- cpm(normdyg,normalized.lib.sizes = TRUE)

# create metadata object
md <- mdcounts %>% select(c("sample_name","genotype","treat","ZT"))
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

md <- treatmap("W")
md <- treatmap("D")

# subset to relevant columns
mspmd <- md %>% select(c(Time,Replicate,Watered,Drought))

# make the design matrix
design <- make.design.matrix(mspmd,degree = 5)
# degree = (number of unique time points) - 1

fit <- p.vector(normcounts, design, Q = 0.05, MT.adjust = "BH")

# find & remove influential genes
tfit<-T.fit(data=fit)

influential<-tfit$influ.info
inf.genenames<-colnames(influential)
nc2<-normcounts[!rownames(normcounts) %in% inf.genenames, ]

# re-run p.vector()
fit2 <- p.vector(nc2,design,Q=0.05,MT.adjust="BH")

# pick k for k-means clustering (I think it will be 7 but let's verify...)
wss<-(nrow(fit2$SELEC)-1)*sum(apply(fit2$SELEC,2,var))

for (i in 2:15){
  wss[i]<- sum(kmeans(fit2$SELEC, centers=i, iter.max=20)$withinss)
} 
plot(1:15, wss, type="b")

pk7kmeans<-see.genes.kh(fit2$SELEC, edesign=fit2$edesign, dis=fit2$dis, cluster.data=1, groups.vector=fit2$groups.vector, cluster.method="kmeans", k=7)

phck7 <- see.genes.kh(fit2$SELEC, edesign=fit2$edesign, dis=fit2$dis, cluster.data=1, groups.vector=fit2$groups.vector, cluster.method="hclust", k=7)
phck5 <- see.genes.kh(fit2$SELEC, edesign=fit2$edesign, dis=fit2$dis, cluster.data=1, groups.vector=fit2$groups.vector, cluster.method="hclust", k=5)
pk5kmeans<-see.genes.kh(fit2$SELEC, edesign=fit2$edesign, dis=fit2$dis, cluster.data=1, groups.vector=fit2$groups.vector, cluster.method="kmeans", k=5)

phck6 <- see.genes.kh(fit2$SELEC, edesign=fit2$edesign, dis=fit2$dis, cluster.data=1, groups.vector=fit2$groups.vector, cluster.method="hclust", k=6)
saveclust(phck6,"../masigpro_results/Yf_hclust_k6_clusters.txt")

saveclust(pk7kmeans,"../masigpro_results/Yf_kmeans_k7_clusters.txt")
saveclust(pk5kmeans,"../masigpro_results/Yf_kmeans_k5_clusters.txt")
saveclust(phck5,"../masigpro_results/Yf_hclust_k5_clusters.txt")
saveclust(phck7,"../masigpro_results/Yf_hclust_k7_clusters.txt")
