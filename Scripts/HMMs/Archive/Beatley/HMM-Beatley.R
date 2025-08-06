#### Beatley HMM ####
rm(list=ls())
source("HMMs/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(tidyverse)

#### Prep Data ####
beatley <- read.csv("Data/Long-Term-Datasets/Beatley/Beatley veg_all_annuals.csv")
beatley <- filter(beatley, YEAR %in% 1963:1967)
select <- c("PLOT", "YEAR", names(which(colSums(beatley[,12:255]) > 0)))
beatley <- beatley[, select]

beatley <- beatley %>% mutate(across(3:71, ~1 * (. > 0))) %>% 
  mutate(across(ACHY:YUSC2, as.numeric)) %>% 
  pivot_longer(names_to = "SpeciesCode", cols = ACHY:YUSC2) %>% 
  pivot_wider(id_cols = c(PLOT, SpeciesCode), values_from = value,  names_from = YEAR)

# Can't use these data, only one plot has successive years