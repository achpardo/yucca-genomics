# Purpose: Run linear modeling on HybridExpress gene numbers to see if there are differences by physiology.
# Author: Anna Pardo
# Date initiated: Apr. 20, 2026

setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/differential_expression/expression_partitioning_downstream/")
library(car)
library(readr)
library(dplyr)
library(ggplot2)

# load data
classcounts <- as.data.frame(read_delim("./hybexp_classcounts.txt",delim = "\t"))

# run overall model for physiology, treat, & Class
omod <- lm(Gene ~ Class + treat + phys + Class*treat + Class*phys + treat*phys +
             Class*treat*phys,data=classcounts)
anova(omod)
# test assumptions: get residuals
res.omod <- resid(omod)
plot(fitted(omod),res.omod,
     main = "HybridExpress full model")
abline(0, 0)
qqnorm(res.omod)
qqline(res.omod)
# run Levene's test for homoscedasticity
result <- leveneTest(Gene ~ interaction(Class,phys,treat),data = classcounts)
print(result)

# try again with log-transformed counts
classcounts$logcount <- log(classcounts$Gene)
# run overall model for physiology, treat, & Class
omod <- lm(logcount ~ Class + treat + phys + Class*treat + Class*phys + treat*phys +
             Class*treat*phys,data=classcounts)
anova(omod)
# test assumptions: get residuals
res.log <- resid(omod)
plot(fitted(omod),res.log,
     main = "HybridExpress full model (log)")
abline(0, 0)
qqnorm(res.log)
qqline(res.log)
# run Levene's test for homoscedasticity
result <- leveneTest(logcount ~ interaction(Class,phys,treat),data = classcounts)
print(result)

# re-run without 3-way interaction
omod2 <- lm(logcount ~ Class + treat + phys + Class*treat + Class*phys + treat*phys,
            data=classcounts)
anova(omod2)
