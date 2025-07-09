#### Height cleaning ####
rm(list = ls())
library(tidyverse)

height <- read.csv("Data/Height/Height_clean/20220805_height_cleaning.csv")
acc <- read.csv("Data/20221216_Seeds_All-Accessions.csv")
species <- read.csv("Data/20211001_Full-Species-List.csv")

height <- filter(height, Choose == "y")

unique(height[!height$Species %in% acc$Species,]$Species)

# get rid of species we aren't using 
cultivated <- unique(acc[acc$site == "GRIN" & acc$type == "cultivated",]$Species)

remove <- unique(height[!height$Species %in% acc$Species,]$Species)

height <- filter(height, !Species %in% cultivated, !Species %in% remove, Species != "Poa secunda")

height$min_max_mean <- (height$max_height_cm + height$min_height_cm)/2

height$Avg_calc_Height_cm <- ifelse(is.na(height$Avg_Height_cm), 
                                    ifelse(is.na(height$min_max_mean),
                                           height$max_height_cm,
                                           height$min_max_mean),
                                    height$Avg_Height_cm)

write.csv(height, "Data/Height/Height_clean/20220805_height_clean.csv", row.names = F)
