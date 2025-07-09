# Merging traits by ID and not by species
rm(list=ls())

### Load Libraries ####
library(plyr)
library(dplyr)
library(tidyverse)
library(ggfortify)
library(gridExtra)
library(Ostats)

#CHECK OUT OSTATS for analysis!

# might want to only choose one variety/subspecies per species
# really should use CP seeds for CP but will need to do that later
# using average for dimorphic seeds

### Enter data ####
acc <- read.csv("Data/20210610_Seeds_All-Accessions.csv")
species <- read.csv("Data/20211001_Full-Species-List.csv")
species <- filter(species, Use != "n")
length(unique(species$Species)) # 293 species total

app <- read.csv("Data/Appendages/20220715_appendages_cleaning.csv")
mass <- read.csv("Data/Seed-mass/Seed-mass_cleaned/20220804_seed-mass_cleaned.csv", na.strings = "")
shape <- read.csv("Data/Shape_Size/20220804_Shape-Size_Cleaned.csv")
disp <- read.csv("Data/Dispersal-vector/Dispersal-vector_cleaned/20220123_dispersal-vector_cleaning.csv") # still in cleaning
set <- read.csv("Data/Settling-time/Settling-time_clean/20220804_Settling-Time_clean.csv")
height <- read.csv("Data/Height/Height_clean/20220805_height_clean.csv")
cn.l <- read.csv("Data/Carbon-nitrogen/Seed-CN.csv")
cn.e <- read.csv("Data/Carbon-nitrogen/MCL-EE_CN.csv")
#oil <- read.csv("Data/Oil/Kew-SID_Oil_Full.csv") # 01-23-2022: this is complete oil but includes protein, search KEW again for protein

### Data Prep ####
#### Appendages ####
app <- filter(app, Type == "primary", Species %in% species$Species, notes != "spikelet")

# for dimorphic seeds i used those with appendages
# presence of appendages is a species level characteristic, I used only CP or other accessions for this variable (if CP had appendages, so do all the other accessions)

unique(app[!app$Species %in% species$Species,]$Species)
unique(species[!species$Species %in% app$Species,]$Species)
length(unique(app$Species)) # 293

#### Mass #### 
mass.sum <- mass %>% filter(chem.morph.updated == "morph")  %>% 
  dplyr::group_by(Species, ID, chem.morph.updated) %>% # there aren't too many duplicates seed wise, with the exception of some from CP and McLaughlin, it would be good to use CP species for the CP data set etc
  dplyr::summarise(mass.mg = mean(Mass.mg, na.rm = T))

unique(species[!species$Species %in% mass.sum$Species,]$Species)
unique(mass.sum[!mass.sum$Species %in% species$Species,]$Species)
length(unique(mass.sum$Species)) # 293

#### Shape & size ####
# Leishman and westoby: shape without appendages
# Thompson et al. 1993: shape with dispersal appendages for some (grasses) but not for others (removed for asters)
# Moles 2000: shape on the persistent part of the diaspore (similar to Thomspon) - this seems like the best method to me
# Funes et al: same here as above, persistent appendages only
#shape <- filter(shape, chem.morph.fun == "morph") # use morphological unit (i.e., with persistent fruit and appendages attached)
shape <- filter(shape, Use == "y")

shape.sum <- ddply(shape, .(Species, ID), summarize, shape = mean(shape), size = mean(size.mm), area = mean(full.area.mm2))

unique(shape.sum[!shape.sum$Species %in% species$Species,]$Species)
unique(species[!species$Species %in% shape.sum$Species,]$Species)
length(unique(shape$Species)) #292

#### Dispersal vector ####
disp <- filter(disp, Quality.control != "n")
none <- "unassisted"
short <- c("ant", "ballistic")
medium <- c("small mammal (internal or caching)", "tumbling", "wind", "bird", "water", "other animal")
long <- c("large mammal (internal)", "adhesive", "human")



disp$ldd <- 0
disp$ldd <- ifelse(disp$main.category %in% short, 1, 
         ifelse(disp$main.category %in% medium, 2, 
                ifelse(disp$main.category %in% long, 3, disp$ldd)))
         
disp <- unique(disp[,c(3:4,10)])
disp <- disp[complete.cases(disp),]
disp$Pres <- 1
disp$cat2 <- disp$main.category
disp$cat2 <- ifelse(disp$main.category == "bird" |
                      disp$main.category == "large mammal (internal)"|
                      disp$main.category == "human", 
                    "ingestion", disp$cat2)

disp$cat2 <- ifelse(disp$main.category == "small mammal (internal or caching)", "caching", disp$cat2)

disp.w <- disp %>% 
  pivot_wider(names_from = "cat2", values_from = "Pres") %>%
  replace(is.na(.), 0) %>%
  dplyr::group_by(Species) %>% 
  dplyr::summarize(across(adhesive:wind, sum)) %>%
  mutate(across(adhesive:wind, as.factor))


disp.sum <- ddply(disp, .(Species), summarize, ldd = ldd[which.max(ldd)], disp = main.category[which.max(ldd)])

length(unique(disp.sum$Species)) # lots are missing :(

remove <- unique(disp.sum[!disp.sum$Species %in% species$Species,]$Species) #these are species i am no longer using

disp.sum <- filter(disp.sum, !Species %in% remove)
disp.w <- filter(disp.w, !Species %in% remove)
disp.sum$ldd2 <- disp.sum$disp
disp.sum$ldd2 <- ifelse(disp.sum$disp == "bird" |
                      disp.sum$disp == "large mammal (internal)"|
                      disp.sum$disp == "human"|
                      disp.sum$disp == "small mammal (internal or caching)", "animal", disp.sum$ldd2)

disp.sum$ldd2 <- ifelse(disp.sum$disp == "tumbling", "wind", disp.sum$ldd2)

disp.sum$ldd2 <- ifelse(disp.sum$disp == "ballistic", "unassisted", disp.sum$ldd2)

disp.sum$ldd2 <- ifelse(disp.sum$disp == "small mammal (internal or caching)", "caching", disp.sum$ldd2)

unique(species[!species$Species %in% disp.sum$Species,]$Species) # 71 species missing, hopefully can get more with genus level data

#### Settling time ####
set <- filter(set, Notes != "floret")
set <- filter(set, Notes != "Floret") # focus on spikelets, the unit that disperses
set.sum <- ddply(set, .(Species, ID), summarize, set.time.mpsec = mean(set.time.mpsec))
unique(set.sum$Species)

unique(set.sum[!set.sum$Species %in% species$Species,]$Species)
unique(species[!species$Species %in% set.sum$Species,]$Species)

rm(set)

#### Height ####
height.sum <- ddply(height, .(Species), summarize, height.cm = mean(Avg_calc_Height_cm)) #293

unique(height.sum[!height.sum$Species %in% species$Species,]$Species)
unique(species[!species$Species %in% height.sum$Species,]$Species) # two species we can't find height for

rm(height)

#### CN ####
cn.l <- merge(cn.l, acc[,c(10,27)], by.x = "Sample.ID", by.y = "ID")
cn  <- rbind(cn.e[,c(2,3,5:19)], cn.l[,c(1,17,2:16)])

cn$prop.C <- cn$C.ug*0.001/cn$sample.wgt.mg
cn$prop.N <- cn$N.ug*0.001/cn$sample.wgt.mg
#cn$prot <- cn$N.ug*0.001
cn$cn <- cn$C.ug/cn$N.ug

cn <- cn %>% 
  dplyr::group_by(Sample.ID, code) %>%
  dplyr::summarise(prop.C = mean(prop.C), prop.N = mean(prop.N), cn = mean(cn))
#cn$prot.prop <- cn$N.ug*0.001/cn$sample.wgt.mg*6.25

remove <- unique(cn[!cn$code %in% species$code,]$code)
cn <- filter(cn, !code %in% remove)

unique(species[!species$code %in% cn$code,]$code) # 4 species left to do: "ANOPEN"  "CHOBREB" "DAUPUS"  "HERHIR" 

cn <- merge(species[,c(3,5)], cn, by = "code")
rm(cn.e, cn.l)

### Merge data ####
full <- merge(mass.sum[,c(1,2,4)], set.sum, by = c("Species", "ID"))
full <- merge(full, height.sum, by = "Species")
full <- merge(full, shape.sum, by = c("Species", "ID"))
full <- merge(full, app[,c(2,4,7,8)], by = "Species")
full <- merge(full, species[,c(3,5)], by = "Species")
full <- merge(full, cn, by.x = c("Species", "ID", "code"), by.y = c("Species", "Sample.ID", "code"), all = T)
full <- merge(full, species[,c(3,7:9)], by = "Species")
full <- merge(full, disp.sum, by = "Species", all.x = T)
full$wing.loading <- full$mass.mg/full$area
  
colnames(full)[10] <- "appendage type"
colnames(full)[12] <- "appendage"
full <- full[,c(1,11,17:19,2:6,8:10,12:16)]

full <- filter(full, Species != "Eschscholzia californica") # not strictly annual
write.csv(full, "Data/20220927_Seed-Traits_ID_cleaning.csv", row.names = F)
