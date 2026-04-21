# Purpose: Run statistics comparing numbers of ASE genes biased each way in different physiological categories and treatments.
# Author: Anna Pardo
# Date initiated: Apr. 17, 2026

setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/differential_expression/")
library(readr)
library(dplyr)
library(ggplot2)

# load counts
asect = as.data.frame(read_delim(file = "./ASE_gttreat_ngenesbiased.txt",delim = "\t"))
asect$ngenes <- asect$`number of gene pairs`

mod <- lm(ngenes ~ ASE,data=asect)
anova(mod)

# put together total number of ASE genes
totals <- asect %>% group_by(CAMphys,genotype,treat) %>% summarize(total_ase_genes = sum(ngenes))

# suite of models to test effect of physiology, genotype, & treatment on total_ase_genes
m1 <- lm(total_ase_genes ~ treat,data=totals)
anova(m1)
m2 <- lm(total_ase_genes ~ CAMphys,data=totals)
anova(m2)

# overall model:
m <- lm(ngenes ~ ASE + CAMphys + genotype + treat +
          ASE*CAMphys + ASE*genotype + ASE*treat + CAMphys*genotype + CAMphys*treat + genotype*treat +
          ASE*CAMphys*genotype + ASE*genotype*treat + CAMphys*genotype*treat +
          ASE*CAMphys*genotype*treat,data=asect)
anova(m)

# genotype is a subset of CAMphys - run with one or the other
cpm <- lm(ngenes ~ ASE + CAMphys + treat + ASE*CAMphys + ASE*treat + CAMphys*treat +
            ASE*CAMphys*treat,data=asect)
anova(cpm)

cpm2 <- lm(ngenes ~ ASE + CAMphys + treat + ASE*CAMphys + ASE*treat + CAMphys*treat,data=asect)
anova(cpm2)
# get residuals
res.cpm2 <- resid(cpm2)
# plot residuals vs. fitted
plot(fitted(cpm2),res.cpm2,
     main = "ASE gene count ~ bias + physiology + treatment + all 2-way interactions")
abline(0, 0) 
# make QQ plot of residuals against the normal distribution
qqnorm(res.cpm2,main = "Normal Q-Q plot of residuals, phys model with 2-way interactions")
qqline(res.cpm2)

# run Levene's test for assumption of homoscedasticity
library(car)
result <- leveneTest(ngenes ~ interaction(ASE,CAMphys,treat),data = asect)
print(result)
