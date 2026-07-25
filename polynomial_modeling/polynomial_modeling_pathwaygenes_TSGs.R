# Purpose: Re-run polynomial modeling with more statistically correct framework (full model only).
# Author: Anna Pardo
# Date initiated: July 22, 2026

library(dplyr)
library(readr)
library(rjson)
library(multcomp)

# set working directory
setwd("//wsl$/Ubuntu/home/leviathan22/Yucca_genomics/rna_insilico_genome/")
# load TPM
tpm <- read_delim("./TPM/Yg_toYgIS_allTPM_correctedmd_over1mil.txt",delim = "\t")

# load physiological genotype categorizations
setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/physiology/")
phys <- fromJSON(file = "./physiological_categories_from_TA.json")
# make a 'phys' column
tpm$phys <- as.character(phys[match(tpm$genotype,names(phys))])
# drop rows with NA for phys
tpm <- tpm %>% filter(phys != "NULL")

# scale TPM before running models
stpm <- as.data.frame(tpm)
row.names(stpm) <- stpm$sample_name
stpm <- stpm %>% select(-c(sample_name,genotype,time,treat,ZT,species,phys))
stpm <- scale(stpm)
stpm <- as.data.frame(stpm)
# add back metadata
stpm$sample_name <- row.names(stpm)
md <- as.data.frame(tpm %>% select(c("sample_name","genotype","phys","treat","ZT")))
stpm_md <- merge(md,stpm)
# drop columns with all NA values
stpm_md <- stpm_md %>% 
  select(
    where(
      ~sum(!is.na(.x)) > 0
    )
  )
stpm_md$phys <- as.factor(stpm_md$phys)
stpm_md$treat <- as.factor(stpm_md$treat)

# load pathway annotations
pathgenes <- as.data.frame(read_csv(file = "../photresp_N_citrate_genes_Yucca.csv"))
clock <- as.data.frame(read_csv("../circadian_light_genes_by_orthology_Yucca.csv"))
cam <- as.data.frame(read_delim("../camgenes_Ya_Yf_orthology_synteny.txt",delim = "\t"))
# put together a list of all the genes in these pathways (vector)
pgenes_vec <- unique(c(clock$GeneID,cam$GeneID,pathgenes$GeneID))
# load 3-way intersection TSGs
tsgs <- read_delim("../masigpro/3way_TSGs_Yg_allphys.txt",delim = "\t",col_names = F)$X1
# create the union of these gene sets
for_analysis <- union(pgenes_vec,tsgs)

########## testing section: work out necessary code ########
gid = for_analysis[1]

# this gene is present in stpm_md
degree = 3

# set up formula
f <- paste(gid,paste0("~ poly(ZT,",as.character(degree),")*treat*phys"))

# run model
fitmod <- lm(f,data = stpm_md)
res <- car::Anova(fitmod,type="III")
resdf <- as.data.frame(car::Anova(fitmod,type="III"))
resdf$Factors <- row.names(resdf)
resdf$GeneID <- gid

# start messing around with glht...
## need to make a matrix showing which comparisons to test
coef(fitmod)
length(coef(fitmod)) # 24 (oy)
summary(glht(fitmod))

# I'm going to use symbolic description option for linfct parameter of the glht() function
## to start with, I just want to know about the physiotype differences from each other
t1 <- glht(fitmod,linfct = c("physCAM = 0",
                             "`physfacultative CAM` = 0"))
summary(t1)
t2 <- glht(fitmod,linfct = mcp(phys = "Tukey"))
summary(t2)
# go with the t2 approach (Tukey-style pairwise comparisons)

s <- summary(t2)
res <- data.frame(
  Comparison = names(s$test$coefficients),
  Estimate   = s$test$coefficients,
  SE         = s$test$sigma,
  t.value    = s$test$tstat,
  p.value    = s$test$pvalues,
  row.names = NULL
)
res

######### Set up functions ########
get_glht_results <- function(glht_obj,
                             test = adjusted("single-step")) {
  s <- summary(glht_obj, test = test)
  
  data.frame(
    Comparison = names(s$test$coefficients),
    Estimate   = s$test$coefficients,
    SE         = s$test$sigma,
    t.value    = s$test$tstat,
    p.value    = s$test$pvalues,
    row.names = NULL
  )
}

run_fullmodel <- function(gid,degree=3,tpmdf = stpm_md){
  if(gid %in% colnames(tpmdf)){
    print(gid)
    print("Running full model...")
    f <- paste(gid,paste0("~ poly(ZT,",as.character(degree),")*treat*phys"))
    fitmod <- lm(f,data = stpm_md)
    resdf <- as.data.frame(car::Anova(fitmod,type="III"))
    resdf$Factors <- row.names(resdf)
    resdf$GeneID <- gid
    
    # run glht
    print("Running post-hoc testing...")
    phres <- glht(fitmod,linfct = mcp(phys = "Tukey"))
    phdf <- get_glht_results(phres)
    phdf$GeneID <- gid
    
    # set up & return a list object
    return(list("model_res" = resdf,"posthoc_res" = phdf))
  }
}

modlist <- list()
phlist <- list()
for(i in 1:length(for_analysis)){
  l <- run_fullmodel(for_analysis[i])
  modlist[[i]] <- l$model_res
  phlist[[i]] <- l$posthoc_res
}
modres <- bind_rows(modlist)
phres <- bind_rows(phlist)

# run FDR for p-values in both dataframes
modres$FDR_p <- p.adjust(modres$`Pr(>F)`,method = "BH")
phres$FDR_p <- p.adjust(phres$p.value,method = "BH")

write.table(modres,file = "../polymod3_modelres_fullmod_TSGs_pathwaygenes.txt",sep = "\t",
            col.names = T,row.names = F)
write.table(phres,file = "../polymod3_posthocres_TSGs_pathwaygenes.txt",sep = "\t",
            col.names = T,row.names = F)
