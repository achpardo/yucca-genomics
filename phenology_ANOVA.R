# Purpose: Run ANOVA on flowering time for different Yucca species.
# Author: Anna Pardo
# Date initiated: June 18, 2026

library(dplyr)
library(readr)

setwd("//wsl$/Ubuntu/home/leviathan22/yucca-genomics/")

# load iNaturalist data
inat <- read_csv(file = "./iNaturalist_obs_all.csv")
# set up a month numbering system
months <- list("Jan"=1,"Feb"=2,"Mar"=3,"Apr"=4,"May"=5,"Jun"=6,"Jul"=7,
               "Aug"=8,"Sep"=9,"Oct"=10,"Nov"=11,"Dec"=12)
# create numerical month column
inat$month_num <- recode(inat$Month,!!!months)

# convert Species to factor
inat$Species <- as.factor(inat$Species)

# subset to just flowering data
flower <- inat %>% filter(Flowering=="Yes")

# try a model
mod <- aov(month_num ~ Species*Latitude, data = flower)
summary(mod)
# Latitude is not significant; drop and re-run as a one-way ANOVA


# run post-hoc test
TukeyHSD(mod)
