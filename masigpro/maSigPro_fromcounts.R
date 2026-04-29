# Purpose: Run maSigPro on the HPC for Y. gloriosa data from each physiological category (categories from titratable acidity data; 3 randomly selected replicates per treatment-ZT combination).
# Author: Anna Pardo
# Date initiated: Apr. 15, 2026

# load modules
library(maSigPro)
library(readr)
library(tools)
library(dplyr)
library(edgeR)

# parse command line args
args = commandArgs(TRUE)
counts = args[1]
outfile = args[2]
# need full paths for BOTH command line arguments

hclust.ap <- function(data,edesign=data$edesign,time.col=1,repl.col=2,group.cols=c(3:ncol(edesign)),names.groups=colnames(edesign)[3:ncol(edesign)],cluster.data=1,groups.vector=data$groups.vector,k=6,distance="cor",agglo.method="ward.D"){
        time = edesign[, time.col]
        repvect = edesign[, repl.col]
        groups = edesign[, group.cols]
        narrays <- length(time)
        clusterdata <- data

        if (!is.null(dim(data))) {
                dat <- as.data.frame(data)
        }

        if (nrow(dat) > 1) {
                dat <- as.data.frame(dat[, (ncol(dat) - length(time) +1):ncol(dat)])
                count.na <- function(x) length(x[is.na(x)])
                NAs <- apply(as.matrix(dat), 1, count.na)
                count.noNa <- function(x) (length(x) - length(x[is.na(x)]))
                dat <- dat[which(apply(as.matrix(dat), 1, count.noNa) >=2), ]
        }
        else {
                NAs <- 1
        }
        kdata <- NULL
        dcorrel<-NULL
        clus2<-NULL
        out <- TRUE
        if (!is.null(clusterdata)){
                #k <- min(k, nrow(dat), na.rm = TRUE)
                dcorrel <- matrix(rep(1, nrow(clusterdata)^2),nrow(clusterdata), nrow(clusterdata)) - cor(t(clusterdata),use = "pairwise.complete.obs")
                clust <- hclust(as.dist(dcorrel), method = agglo.method)
                c.algo.used = paste("hclust", "cor",agglo.method, sep = "_")
                cut <- cutree(clust,k=k)
                groups <- as.matrix(groups)
                colnames(groups) <- names.groups
        }
        OUTPUT <- list(cut, c.algo.used, groups, dcorrel)
	names(OUTPUT) <- c("cut", "cluster.algorithm.used", "groups", "correl")
        OUTPUT}

# data normalization function
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

# function to construct watered & drought columns for metadata
treatmap <- function(t,md){
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


# metadata setup function
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
	md <- treatmap("W",md)
	md <- treatmap("D",md)

	# subset to relevant columns
	mspmd <- md %>% select(c(Time,Replicate,Watered,Drought))
	# make the design matrix
	design <- make.design.matrix(mspmd,degree = 5)
	return(design)
}

# function to run maSigPro
run_masigpro <- function(normcounts,design,rdfname){
	fit <- p.vector(normcounts, design, Q = 0.05, MT.adjust = "BH")
	# find & remove influential genes
	tfit<-T.fit(data=fit)
	influential<-tfit$influ.info
	inf.genenames<-colnames(influential)
	nc2<-normcounts[!rownames(normcounts) %in% inf.genenames, ]
	# save an RData object (to save waiting time on local)
	save.image(rdfname)
	# re-run p.vector()
	fit2 <- p.vector(nc2,design,Q=0.05,MT.adjust="BH")
	return(fit2)
}

run_all <- function(counts,rdfname){
	tmm <- make_tmm(counts)
	des <- metadata_setup(counts)
	fit2 <- run_masigpro(tmm,des,rdfname)
	prof <- hclust.ap(fit2$SELEC, edesign=fit2$edesign, dis=fit2$dis, cluster.data=1,groups.vector=fit2$groups.vector, k=6)
	# extract genes in clusters
	clust <- as.data.frame(prof$cut)
	clust$GeneID <- row.names(clust)
	clust$cluster <- clust$`prof$cut`
	clust <- clust %>% select(c(cluster,GeneID))
	return(clust)
}

rdfname <- paste0(counts,"_for_troubleshooting.RData")
# load counts data
counts <- as.data.frame(read_delim(counts,delim="\t"))

# run maSigPro time-structured gene ID
clust <- run_all(counts)

# save output file
write.table(clust,outfile,sep="\t",col.names=T,row.names=F)
