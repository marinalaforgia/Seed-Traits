# Merging traits
rm(list=ls())

### Load Libraries ####
library(plyr)
library(dplyr)
library(tidyverse)
library(ggfortify)
library(gridExtra)

# using average for dimorphic seeds

### Enter data ####
acc <- read.csv("Data/20230530_Seeds_All-Accessions.csv")
length(unique(acc$Species))

species <- read.csv("Data/20211001_Full-Species-List.csv")
species <- filter(species, Use != "n")
length(unique(species$Species)) # 293 species total
length(unique(species$Species.no.ssp)) #278 total without subspecies

acc <- acc[acc$Species %in% species$Species,]
length(unique(acc$Species))

## Individual traits ##
app <- read.csv("Data/Appendages/20220715_appendages_cleaning.csv")
mass <- read.csv("Data/Seed-mass/Seed-mass_cleaned/20230323_seed-mass_cleaned.csv", na.strings = "")
shape <- read.csv("Data/Shape_Size/20230323_Shape-Size_Cleaned_hand-curated.csv")
disp <- read.csv("Data/Dispersal-vector/Dispersal-vector_cleaned/20230522_dispersal-vector_clean.csv")
set <- read.csv("Data/Settling-time/Settling-time_clean/20221216_Settling-Time_clean.csv")
height <- read.csv("Data/Height/Height_clean/20220805_height_clean.csv")
cn.l <- read.csv("Data/Carbon-nitrogen/Seed-CN.csv")
cn.e <- read.csv("Data/Carbon-nitrogen/MCL-EE_CN.csv")
oil <- read.csv("Data/Oil/Kew-SID_Oil_Full.csv") # 01-23-2022: this is complete oil but includes protein, search KEW again for protein
int <- read.csv("Data/Internal-traits/20230522_Internal-Morphology_cleaned.csv")
coat <- read.csv("Data/Coat-permeability/20230324_Seed-Coat-Permeability_cleaned.csv")

### Data Prep ####
#### Appendages ####
app <- filter(app, Type == "primary", Species %in% species$Species, notes != "spikelet")

app[app$Species == "Malacothrix fendleri",]$descriptor <- "pappus" #eflora indicates this should have a pappus, likely removed prior to obtaining

app[app$Species == "Palafoxia arida var. arida",]$descriptor <- "pappus"


# for dimorphic seeds i used those with appendages
# presence of appendages is a species level characteristic, I used only CP or other accessions for this variable (if CP had appendages, so do all the other accessions)

unique(app[!app$Species %in% species$Species,]$Species)
unique(species[!species$Species %in% app$Species,]$Species)
length(unique(app$Species)) # 293

#### Mass #### 
mass.sum <- mass %>% filter(chem.morph.updated == "morph")  %>% 
  dplyr::group_by(Species, site, chem.morph.updated, ID) %>% # there aren't too many duplicates seed wise, with the exception of some from CP and McLaughlin, it would be good to use CP species for the CP data set etc
  dplyr::summarise(morph.mass.mg = mean(Mass.mg, na.rm = T))

mass.sum.chem <- mass %>% filter(chem.morph.updated == "chem")  %>% 
  dplyr::group_by(Species, site, chem.morph.updated, ID) %>% # there aren't too many duplicates seed wise, with the exception of some from CP and McLaughlin, it would be good to use CP species for the CP data set etc
  dplyr::summarise(mass.mg = mean(Mass.mg, na.rm = T))

colnames(mass.sum.chem)[5] <- "chem.mass.mg"

mass.sum <- merge(mass.sum[,-3], mass.sum.chem[,-3], by = c("Species", "site", "ID"), all = T)

mass.sum$chem.mass.mg <- ifelse(is.na(mass.sum$chem.mass.mg), mass.sum$morph.mass.mg, mass.sum$chem.mass.mg)

unique(species[!species$Species %in% mass.sum$Species,]$Species)
unique(mass.sum[!mass.sum$Species %in% species$Species,]$Species)
length(unique(mass.sum$Species)) # 293

rm(mass, mass.sum.chem)

#### Shape & size ####
# Update 12-16-2022: went through by hand to choose which ones we want to keep, decided on persistent appendages only
# Leishman and westoby: shape without appendages
# Thompson et al. 1993: shape with dispersal appendages for some (grasses) but not for others (removed for asters)
# Moles 2000: shape on the persistent part of the diaspore (similar to Thomspon) - this seems like the best method to me
# Funes et al: same here as above, persistent appendages only

shape.sum <- filter(shape[,-c(6,7)], Use == "y") # use morphological unit (i.e., with persistent fruit and appendages attached; this removes appendages for asters and other appendages that readily fall off)

shape.sum <- merge(shape.sum, acc[,c(1,7,25)], by = c("Species", "ID"), all.x = T, all.y = F)

shape.sum <- ddply(shape.sum, .(Species, site, ID), summarize, shape = mean(shape), size.mm = mean(size.mm), width = mean(width.mm))

unique(shape.sum$Species) #293

# wing loading requires area that includes pappus, use morph
wing <- filter(shape, chem.morph.fun == "morph")
wing <- merge(wing, mass.sum, by = c("Species", "ID"), all = T)
wing$wing.loading <- wing$morph.mass.mg/wing$full.area.mm2
unique(wing$Species) #293

wing <- filter(wing, !(ID == "CBG-20891" & Use == "n")) # this rep did not have the app despite being labeled as "morph"

shape.sum[!shape.sum$Species %in% wing$Species,]$Species #3 missing due to mismatch in which acc we did mass on and which accession we did area on MEDPOL, PHAACU, THYCUR

# MEDPOL - area and mass on two diff accessions, but both still from MCL
wing[wing$Species == "Medicago polymorpha",]$wing.loading  <- 165.51/158.1747

# THYCUR - area and mass on two diff accessions, one CBG (area) other MCL (mass)
wing[wing$Species == "Thysanocarpus curvipes",]$wing.loading  <- 13.59/49.91892

# PHAACUT - all from GRIN
wing[wing$Species == "Phaseolus acutifolius var. tenuifolius",]$wing.loading  <- ((220.28 + 207.82)/2)/31.76495

wing <- wing[,c(1,2,15)]
wing <- wing[complete.cases(wing),]

length(unique(wing$Species))

#explore both chem and morph shapes
# shape.c <- filter(shape, chem.morph.fun == "chem")
# shape.m <- filter(shape, chem.morph.fun == "morph")
# 
# colnames(shape.c)[6:8] <- c("size.mm.chem", "shape.chem", "full.area.mm2.chem" )
# colnames(shape.m)[6:8] <- c("size.mm.morph", "shape.morph", "full.area.mm2.morph" )
# 
# shape2 <- rbind(shape.c, shape.m)
# 
# shape.sum <- ddply(shape2, .(Species), summarize, shape = mean(shape), size = mean(size.mm), area = mean(full.area.mm2))
# 
# 
# shape.c.sum <- ddply(shape.c, .(Species), summarize, shape.c = mean(shape.chem), size.mm.c = mean(size.mm.chem), area.mm2.c = mean(full.area.mm2.chem))
# 
# shape.m.sum <- ddply(shape.m, .(Species), summarize, shape.m = mean(shape.morph), size.mm.m = mean(size.mm.morph), area.mm2.m = mean(full.area.mm2.morph))
# 
# shape.sum <- merge(shape.m.sum, shape.c.sum, by = "Species", all = T)
# 
# shape.sum$shape.c <- ifelse(is.na(shape.sum$shape.c), shape.sum$shape.m, shape.sum$shape.c)
# shape.sum$size.mm.c <- ifelse(is.na(shape.sum$size.mm.c), shape.sum$size.mm.m, shape.sum$size.mm.c)
# shape.sum$area.mm2.c <- ifelse(is.na(shape.sum$area.mm2.c), shape.sum$area.mm2.m, shape.sum$area.mm2.c)
# 
# unique(shape.sum[!shape.sum$Species %in% species$Species,]$Species)
# unique(species[!species$Species %in% shape.sum$Species,]$Species)
# length(unique(shape$Species)) #293
# 
# rm(shape, shape.c, shape.m, shape.m.sum, shape.c.sum)


#### Settling time ####
set <- filter(set, Notes != "floret")
set <- filter(set, Notes != "Floret") # focus on spikelets, the unit that disperses
set <- merge(set, acc[,c(1,7,25)], by = c("Species", "ID"), all.x = T, all.y = F)

set <- filter(set, ID != "CBG-22362") # removing chem unit NAVPUB as dispersal is at the morph unit 

set.sum <- ddply(set, .(Species, site, ID), summarize, set.time.mpsec = mean(set.time.mpsec))
unique(set.sum$Species)

unique(set.sum[!set.sum$Species %in% species$Species,]$Species)
unique(species[!species$Species %in% set.sum$Species,]$Species)

rm(set)

#### Height ####
height.sum <- ddply(height, .(Species), summarize, height.cm = mean(Avg_calc_Height_cm)) #293

unique(height.sum[!height.sum$Species %in% species$Species,]$Species)
unique(species[!species$Species %in% height.sum$Species,]$Species) 

rm(height)

#### Dispersal vector ####
# dispersal distance classes based on Lasosova et al 2023
# class 1: unassisted and under 30 cm
# class 2: unassisted and over 30 cm
# class 3: ants
# class 4: wind (tumbleweed)
# class 5: wind (pappus)
# class 6: animal (includes adhesion and ingestion)
# class 7: human 

disp <- filter(disp, Quality.control != "n", Quality.control != "", main.category != "water", Quality.control.notes != "not primary") # removing water because this is associated with washes (which may be important but arent primary vectors)

disp <- merge(disp, height.sum, by = "Species")

unique(disp$main.category)

animal <- c("small mammal (internal or caching)", "bird", "other animal", "large mammal (internal)")

disp$main.category <- recode_factor(disp$main.category, `small mammal (internal or caching)` = "animal", bird = "animal", `other animal` = "animal", `large mammal (internal)` = "animal")

unknown <- "unknown"
class.1.2 <- c("unassisted", "ballistic")
class.3 <- c("ant", "other insect")
class.4 <- "tumbling"
class.5 <- "wind"
class.6 <- c("animal", "adhesive")
class.7 <- "human"


disp$ldd <- NA
disp$ldd <- ifelse(disp$main.category %in% class.1.2 & disp$height.cm <= 30, 1, 2)
disp$ldd <- ifelse(disp$main.category %in% class.3, 3, disp$ldd)
disp$ldd <- ifelse(disp$main.category %in% class.4, 4, disp$ldd)
disp$ldd <- ifelse(disp$main.category %in% class.5, 5, disp$ldd)      
disp$ldd <- ifelse(disp$main.category %in% class.6, 6, disp$ldd)  
disp$ldd <- ifelse(disp$main.category %in% class.7, 7, disp$ldd)  

disp <- unique(disp[,c(1,4,12)])

disp.sum <- ddply(disp, .(Species), summarize, ldd.all = ldd[which.max(ldd)], disp.cat.all = main.category[which.max(ldd)])
disp.sum2 <- ddply(disp[disp$ldd != 7, ], .(Species), summarize, ldd.natural = ldd[which.max(ldd)], disp.cat.nat = main.category[which.max(ldd)])

disp <- merge(disp.sum, disp.sum2, by = "Species")

remove <- unique(disp.sum[!disp.sum$Species %in% species$Species,]$Species) #these are species i am no longer using

disp.sum <- filter(disp.sum, !Species %in% remove)

unique(species[!species$Species %in% disp.sum$Species,]$Species) 

#### CN ####
cn.l <- merge(cn.l, acc[,c(1,25,9)], by.x = "Sample.ID", by.y = "ID")
cn.e$site <- "MCL-EE"
cn.e <- merge(acc[,c(1,25,9)], cn.e, by = c("code", "site"), all.x = F, all.y = T)

colnames(cn.l)[1] <- "ID"

cn  <- rbind(cn.e[,c(3,2,1,7:21)], cn.l[,c(1,17,18,2:16)])

cn$prop.C <- cn$C.ug*0.001/cn$sample.wgt.mg
cn$prop.N <- cn$N.ug*0.001/cn$sample.wgt.mg
cn$cn <- cn$C.ug/cn$N.ug

cn <- cn %>% 
  dplyr::group_by(code, site, ID) %>%
  dplyr::summarise(prop.C = mean(prop.C), prop.N = mean(prop.N), cn = mean(cn))
#cn$prot.prop <- cn$N.ug*0.001/cn$sample.wgt.mg*6.25

remove <- unique(cn[!cn$code %in% species$code,]$code)
cn <- filter(cn, !code %in% remove)

unique(species[!species$code %in% cn$code,]$code) # 1 species we did not have enough seed - "HERHIR" 

cn <- merge(species[,c(2,5)], cn, by = "code")

# AMSMEN has a suspiciously low prop.C compared to other Amsinckia species, replace with mean
cn[cn$Species == "Amsinckia menziesii",]$prop.C <- mean(cn$prop.C, na.rm = T)

rm(cn.e, cn.l)

#### Internal traits ####
# int$coat.all <- ifelse(!is.na(int$thickness.coat.mm), 
#                        int$thickness.coat.mm, 
#                        ifelse(!is.na(int$thickness.fruit.mm), 
#                               int$thickness.fruit.mm,
#                               int$thickness.both.mm))

int[is.na(int$thickness.coat.mm),]$thickness.coat.mm <- 0 
int[is.na(int$thickness.fruit.mm),]$thickness.fruit.mm <- 0 
int$thickness.both.mm <- ifelse(is.na(int$thickness.both.mm), int$thickness.coat.mm + int$thickness.fruit.mm, int$thickness.both.mm)

int$E.S <- ifelse(!is.na(int$segmented.embryo.length), int$segmented.embryo.length/int$seed.length.mm, int$embryo.length.mm/int$seed.length.mm) 

int <- merge(int, acc[,c(1,7,25)], by = c("Species", "ID"), all.x = T, all.y = F)

int.sum <- int %>%
  dplyr::group_by(Species, site, ID) %>%
  dplyr::summarise(E.S = mean(E.S, na.rm = T), 
                   coat.thick = mean(thickness.coat.mm, na.rm = T),
                   fruit.thick = mean(thickness.fruit.mm, na.rm = T),
                   both.thick = mean(thickness.both.mm, na.rm = T),
                   mucilage = mean(mucilage, na.rm = T),
                   fleshy.end = mean(fleshy.endosperm, na.rm = T),
                   starchy.end = mean(starchy.endosperm, na.rm = T))


#### Seed Coat permeability ####
coat$mass.diff.mg <- coat$end.mass.mg - coat$start.mass.mg
coat$perc.diff <- coat$mass.diff.mg/coat$start.mass.mg * 100

coat <- merge(coat, acc[,c(1,7,25)], by = c("Species", "ID"), all.x = T, all.y = F)

coat.sum <- coat %>%
  #filter(perc.diff >= 0) %>% # only a couple that lost weight
  dplyr::group_by(Species, site, ID) %>%
  dplyr::summarise(coat.perm.perc = mean(perc.diff, na.rm = T))


### Merge data ####
full <- merge(mass.sum, set.sum, by = c("Species", "site", "ID"), all = T)
full <- merge(full, cn[,2:6], by = c("Species", "site", "ID"), all = T)
full <- merge(full, height.sum, by = "Species", all.x = T, all.y = F)
full <- merge(full, shape.sum, by = c("Species", "site", "ID"), all = T)
full <- merge(full, coat.sum, by = c("Species", "site", "ID"), all = T)
full <- merge(full, int.sum, by = c("Species", "site", "ID"), all = T)
full <- merge(full, app[,c(2,4,7,8)], by = "Species")
full <- merge(full, disp, by = "Species", all.x = T)
full$coat.thick.per.width <- full$both.thick/full$width
full <- merge(full, wing[,c(2,3)], by = "ID", all.x = T)
colnames(full)[21] <- "appendage type"
colnames(full)[23] <- "appendage"

full <- merge(full, species[,c(2,4,5,6, 8:10)], by = "Species")
full <- merge(full, acc[,c(1,30)], by = "ID", all.x = T, all.y = F) # add in AI data


full <- full[,c(2,1,3, 30:35, 21, 23, 25, 27, 4:17, 28, 29, 24, 26, 18:20, 22, 36)]

# oil 
colnames(oil)[3:5] <- c("oil.prop.kew", "prot.prop.kew", "carb.prop.kew")
full <- merge(full, oil[,c(1,3:5)], by.x = "Species", by.y = "species", all.x = T, all.Y = F)

full[full$Species == "Micropus californicus" & full$site == "SFREC",]$chem.mass.mg <- NA # this number is not the actual chem unit, that is only in the CPB accession
full[full$Species == "Micropus californicus" & full$site == "CBG",]$set.time.mpsec <- NA # similarly the unit that disperses is the morph unit so settling time of chem unit meaningless
full[full$Species == "Aegilops triuncialis" & full$ID == "SFREC-20",]$chem.mass.mg <- NA # this unit is the spikelet, thus does not have a chem mass

full[full$Species == "Chorizanthe fimbriata" & full$ID == "CBG-21158",]$chem.mass.mg <- NA

full[full$Species == "Chorizanthe staticoides" & full$ID == "SS-11",]$chem.mass.mg <- NA

# arizonicum and neomexicanum synonyms
full[full$new.code == "CHEARI",]$new.code <- "CHENEO"

write.csv(full, "Data/20230801_Seed-Traits_clean_ID.csv", row.names = F)

### Summarize traits ####
# by species (removing subspecies) and site
full <- read.csv("Data/20230801_Seed-Traits_clean_ID.csv") # edited keep column to indicate which accessions should be combined/averaged and which to keep separate (like carrizo)

# group by species, site, and whether or not to keep separate or combine, then calculate average seed trait for that site (averages across instances where chemical traits were done on one acc due to lack of substantial seeds but every other traits were done on a diff accession, keeps CP separate from MCL, except for a couple instances where C & N were done on other accessions)

full.comb <- full %>% 
  filter(keep != "discard") %>%
  group_by(Species.no.ssp, keep, site) %>%
  dplyr::summarize(across(morph.mass.mg:carb.prop.kew, ~ mean(.x, na.rm = TRUE)))

full.comb$coat.thick.per.width <- ifelse(is.nan(full.comb$coat.thick.per.width), full.comb$both.thick/full.comb$width, full.comb$coat.thick.per.width)

full.unique <- unique(full[,c(6,8:15)])

#full.unique[duplicated(full.unique$Species.no.ssp),]

full <- merge(full.unique, full.comb[,-c(2,12)], by = "Species.no.ssp", all = T)

colnames(full)[1] <- "Species"

### Combine with AI range data ####
# herb.AI <- readRDS("Data/Herbarium-Occurences/Herbarium_occ_AI_sum.RDS")
# 
# full <- merge(full, herb.AI, by.x = "Species", by.y = "Marina.species", all = T)

write.csv(full, "Data/20230801_Seed-Traits_clean_site.csv", row.names = F)

####
ggplot(traits, aes(x = prop.C, y = carb.prop)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlim(0.375,0.575) # why in the world are proportion carbon as measured at SIF and proportion carbohydrate in the KEW database NEGATIVELY correlated?

ggplot(traits, aes(x = prop.N*6.25, y = prot.prop)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) + 
  geom_smooth(method = "lm") # this makes sense

ggplot(traits, aes(x = prop.C, y = oil.prop)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1) + 
  xlim(0.375,0.575) # carbon is a building block of lipids and carbs, i guess in this dataset, carbon is a better predictor of oil than of carbs

#carbon is a building block in both lipids and carbohydrates

# which species are exclusively arid, which are exclusively semi-arid and which are in both environments?

# use cch2 to answer that?

# download cch2 for each species