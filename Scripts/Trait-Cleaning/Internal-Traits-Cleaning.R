# Internal traits cleaning ##
library(tidyverse)
library(stringr)

rm(list=ls())

int1 <- read.csv("Data/Internal-traits/20230330_Internal-Morphology_entered.csv")
int2 <- read.csv("Data/Coat-permeability/20230324_Seed-Coat-Permeability_entered.csv")

int1 <- int1[,c(1:5, 7, 8, 9, 11, 12, 13, 17:24)]
  

species <- read.csv("Data/20211001_Full-Species-List.csv")
species <- filter(species, Use != "n")

# internal traits
int1$Species <- recode_factor(int1$Species, 'Melilotus indica' = "Melilotus indicus")
unique(int1[!int1$Species %in% species$Species,]$Species) # no mispellings

str(int1)
int1$mucilage <- ifelse(int1$mucilage == "y", 1, 0)
int1$endosperm.present <- ifelse(int1$endosperm.present == "y", 1, 0)
int1$starchy.endosperm <- ifelse(int1$starchy.endosperm == "y", 1,
                        ifelse(int1$starchy.endosperm == "n", 0, NA))
int1$fleshy.endosperm <- ifelse(int1$fleshy.endosperm == "y", 1,
                        ifelse(int1$fleshy.endosperm == "n", 0, NA))

str(int1)

int1 <- int1[,-c(11,15,19)]

int1$thickness.both.mm <- ifelse(!is.na(int1$thickness.coat.mm) & !is.na(int1$thickness.fruit.mm), int1$thickness.coat.mm + int1$thickness.fruit.mm, int1$thickness.both.mm)

int1$embryo.length.mm <- ifelse(!is.na(int1$segmented.embryo.length), int1$seed.length.mm, int1$embryo.length.mm)

write.csv(int1, "Data/Internal-traits/20230522_Internal-Morphology_cleaned.csv", row.names = F)

pics <- list.files("Data/Internal-traits/Photos/")
pics <- str_extract(pics, "[^_]+")

unique(pics) #298
species[!species$code %in% pics,]$code # none!
unique(pics[!pics %in% species$code]) # all no longer in use

test <- table(pics)

# seed perm
int2 <- int2[,c(2:5,11:12)] # went through notes
int2 <- filter(int2, !is.na(end.mass.mg))

write.csv(int2, "Data/Coat-permeability/20230324_Seed-Coat-Permeability_cleaned.csv", row.names = F)

## found some misentered thicknesses
int1 <- read.csv("Data/Internal-traits/20230522_Internal-Morphology_cleaned.csv")

tmp1 <- filter(int1, !is.na(thickness.both.mm))
ggplot(tmp1[401:430,], aes(x = ID, y = thickness.both.mm)) +
  geom_boxplot()
