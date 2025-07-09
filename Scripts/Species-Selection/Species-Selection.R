#### Seed selection for trait studies ###
rm(list=ls())

#### Libraries ####
library(plyr)
library(tidyverse)

#### Species databases ####
## USDA Plants California
CA.usda <- read.csv("Projects/Seed-Traits/Datasets/USDA-Plants_CA_Annuals.csv") # dataset contains some alternate names

## 1. Cal-IPC: California Invaders
cal.ipc <- read.csv("Projects/Seed-Traits/Datasets/Cal-IPC_Annual.csv") #downloaded from CalFlora, filter options: annual, cal-ipc invaders, sorted by family; note there are many non-natives not included here
cal.ipc$Species <- as.character(cal.ipc$Species)
cal.ipc <- filter(cal.ipc, collect. == "Y")
cal.inv <- read.csv("Projects/Seed-Traits/Datasets/cal-ipc_inventory.csv")
cal.inv <- filter(cal.inv, Scientific.name %in% cal.ipc$Species, )

## 2. Calflora: annuals in CA downloaded form calflora, filter options: annual, native, sorted by family
calflora <- read.csv("Projects/Seed-Traits/Datasets/Calflora_nat-ann.csv")
calflora$Species <- as.character(calflora$Species)
calflora$Family <- as.character(calflora$Family)
CA <- rbind(cal.ipc, calflora)
CA <- merge(CA, CA.usda, by = c("Species", "Family", "group", "nat.inv"), all.x = T)
CA <- CA[,c(5,1:4)]
rm(cal.ipc, calflora)

## 3. JR: Jasper Ridge Seeds
JR <- read.csv("Projects/Seed-Traits/Datasets/JR_annual_species.csv")

## 4. SD: Sonoran Desert Seeds
SD <- read.csv("Projects/Seed-Traits/Datasets/Sonoran-Desert_Species-List_2016.csv")

## 5. Arizona Annuals: USDA plants database, family, annual, invasive status
AZ <- read.csv("Projects/Seed-Traits/Datasets/Arizonia_Annuals.csv")

## 6. Kew Seed Database for Oil content
oil <- read.csv("Projects/Seed-Traits/Datasets/Kew-SID_Oil_Full.csv")
oil$Species <- paste(oil$Genus, oil$Species, oil$var.ssp, oil$var.ssp.name)
oil$Species <- trimws(oil$Species)
oil <- oil[,c(2,5,6,7)]
oil <- filter(oil, Species %in% CA$Species | Species %in% AZ$Species | Species %in% CA.usda$Species) #378
oil.ca <- merge(oil, CA, by = "Species", all.x = F)
oil.az <- merge(oil, AZ, by = "Species", all.x = F)
oil.usda <- merge(oil, CA.usda, by = "Species", all.x = F)
oil <- rbind(oil.ca, oil.az)
oil <- rbind(oil, oil.usda)
oil <- distinct(oil)
oil$kew.oil <- "y"
rm(oil.ca, oil.az, oil.usda)

# working species list
pos.spp <- read.csv("Projects/Seed-traits/Datasets/20210112_Possible-Species.csv")
#oil.pos <- filter(oil, Species %in% pos.spp$Species)

#### Seed Acquisition ####
#Notes : DBG's collected is over 20 years old, often from only one plant (so very few seeds) and only 30 species, not worth the trouble; Kew Millenium seed bank only gives out 60 seeds per accession which isn't much; Daniel's collections of invaders could supplement mine, though really his most useful species is brassica tournefortii

## CBG - California botanic garden FKA Rancho Santa Ana Botanic Garden
cbg.stat <- read.csv("Projects/Seed-Traits/Datasets/RSABG-Online-Seed-List-Jan-2020.csv") # CBG file with collection status (filters out seeds not available for request: N/A, SR/ST/SE, FE/FT)
cbg <- read.csv("Projects/Seed-Traits/Datasets/RSABG-Dryland-annuals_2021.csv") # file has been clipped to dryland ecoregions and filtered to annuals
cbg <- merge(cbg, CA, by = "Species", all.x = T)

## S & S: use this for oil, embryo morphology
s.s <- read.csv("Projects/Seed-Traits/Datasets/S-and-S_Seeds.csv")
s.s$Site <- "ss"
#filter(s.s, Species %in% oil$Species)
s.s$Site <- ifelse(s.s$Species %in% cbg$Species, "both", "ss")
cbg$Site <- ifelse(cbg$Species %in% s.s$Species, "both", "cbg")
cbg <- merge(cbg, s.s, by = c("Species", "Family", "group", "nat.inv", "Site"), all.x = T, all.y = T)
cbg <- merge(cbg, oil, by = c("Species", "Family", "group", "nat.inv"), all.x = T, all.y = T)
cbg[is.na(cbg$kew.oil), ]$kew.oil <- "n"

## GRIN: Arizona species
grin <- read.csv("Projects/Seed-Traits/Datasets/USDA-GRIN-AZ.csv")
grin$Received <- as.numeric(as.character(grin$Received))
grin <- filter(grin, !Species %in% test$Species) # species that dont overlap with CBG
grin <- filter(grin, Species %in% AZ$Species | Species %in% CA$Species | Species %in% CA.usda$Species, Received > 1995)
length(grin.sp <- unique(grin$Species)) #57

write.csv(grin, "Projects/Seed-Traits/Datasets/USDA-GRIN-AZ-annuals.csv", row.names = F)

grin.sp <- merge(data.frame(Species = grin.sp), AZ, by = "Species")
grin.sp$Site <- "grin"
grin.sp <- merge(grin.sp[,-2], oil[,c(1,9)], by = "Species", all.x = T)
grin.sp[is.na(grin.sp$kew.oil),]$kew.oil <- "n"

## McLaughlin (Elise)
McL <- read.csv("Projects/Seed-Traits/Datasets/McLaughlin_Species_EE.csv")
McL <- merge(McL, CA.usda, by = "Species", all.x = T)
write.csv(McL, "Projects/Seed-Traits/Datasets/Elise_Species.csv", row.names = F)

## Daniel Winkler (Desert invaders)

#### CBG ####
# Cannot request more than 10% of an accession
cbg$seed_.15g <- cbg$seed.per.g*0.15
#cbg.sum <- ddply(cbg[cbg$date.store > 1995,], .(Species), summarize, tot.seed = sum(seed.num), tot.wt = sum(tot.wt), seed.per.g = mean(seed.per.g))
cbg$seed.num.1 <- 0.1*cbg$seed.num
cbg2 <- filter(cbg, seed.num.1 > 300, date.store > 1995)
cbg2$true <- ifelse(cbg2$seed_.15g < cbg2$seed.num.1, T, F)
#cbg2 <- filter(cbg2, true == T)
write.csv(cbg2, "Projects/Seed-Traits/Datasets/Orders/CBG_order_LaForgia.csv", row.names = F)

test <- filter(cbg, 0.1*seed.num > seed_.2g, date.store > 1995)
length(unique(cbg2$Species)) #180 species; if they allow 20% of an accession i can get it to 229 species


#180 species from CBG
length(unique(test$Species))

cbg.sum$seed_.2g <- cbg.sum$seed.per.g*0.2
cbg.sum <- filter(cbg.sum, 0.1*tot.seed > seed_.2g)
cbg2 <- filter(cbg[,c(1:5,11,23)], Species %in% cbg.sum$Species, date.store > 1995) # what are storage conditions like?
cbg2 <- distinct(cbg2[,c(1:5,7)])

#20 species at McLaughlin from Elise that I can supplement with SS for other traits
McL <- filter(McL, collected == "y")
McL.ss <- filter(McL, Species %in% s.s$Species)
McL.ss$Site <- "McL.ss"
McL.ss <- merge(McL.ss[,c(1,6,7,8,9)], oil[,c(1,9)], by = "Species", all.x = T)
McL.ss[is.na(McL.ss$kew.oil),]$kew.oil <- "n"
McL.only <- filter(McL.cbg, !Species %in% cbg2$Species)
McL.only$Site <- "McL.E"
write.csv(McL.only,"Projects/Seed-Traits/McLaughlin_Elise-subset.csv", row.names = F)

# Grin has 57 species

# Daniel has 10 species

full <- rbind(cbg2, grin.sp)
full <- rbind(full, McL.ss)
write.csv(full, "Projects/Seed-Traits/20210203_Possible-Species.csv", row.names = F)
pos.spp <- merge(pos.spp, CA.usda[,2:5], by = c("Species", "group", "nat.inv"), all.x = T)
pos.spp2 <- filter(pos.spp, !Species %in% full$Species)
pos.spp <- pos.spp[,c(1,6,2,3)]
pos.spp <- merge(pos.spp, oil[,c(1,9)], by = "Species", all.x = T)
pos.spp[is.na(pos.spp$kew.oil),]$kew.oil <- "n"

#####
full <- read.csv("Projects/Seed-Traits/20210203_Possible-Species.csv")
# full <- merge(full, pos.spp, by = c("Species", "Family", "group", "nat.inv", "kew.oil"), all = T)
write.csv(full, "Projects/Seed-Traits/20210203_Possible-Species.csv", row.names = F)

SD <- filter(SD, !Species %in% full$Species)
SD <- filter(SD, !Alt.name %in% full$Species)
filter(SD, Species %in% s.s$Species)
JR <- filter(JR, !Species %in% full$Species)
filter(JR, Species %in% s.s$Species)
ss.no <- filter(s.s, !Species %in% full$Species)


kew.mass <- read.csv("Projects/Seed-Traits/Datasets/Archive/Kew_SID_Seed-Mass_Full.csv")
kew.mass$Species <- paste(kew.mass$Genus, kew.mass$Species, kew.mass$var.ssp, kew.mass$var.ssp.id)
kew.mass <- kew.mass[,c(2,5)]
kew.mass$Species <- trimws(kew.mass$Species)
kew.mass1 <- filter(kew.mass, Species %in% full$Species)
kew.massCA <- filter(kew.mass, Species %in% CA$Species)
kew.massAZ <- filter(kew.mass, Species %in% AZ$Species)
kew.mass2 <- filter(full, !Species %in% kew.mass$Species)

# Leptosiphon montanus, Hemizonia fasciculata, Erythranthe guttata, Festuca microstachys, Festuca myuros are missing due to name change

full <- merge(full, kew.mass, by = "Species", all.x = T)

#### 
rm(list=ls())
cbg <- read.csv("Projects/Seed-Traits/Datasets/RSABG-Dryland-annuals_2021.csv") 
full <- read.csv("Projects/Seed-Traits/Datasets/Orders/cbg_order.csv")

cbg2 <- filter(cbg, date.store > 1995)

full2 <- filter(full, !full$Species %in% cbg2$Species)

cbg2$seed_10pct <- cbg2$seed.num*0.1
cbg2 <- filter(cbg2, seed_10pct > 300)
cbg2$seed.2 <- cbg2$seed.per.g*0.15
cbg2$true <- ifelse(cbg2$seed_10pct > cbg2$seed.2, T, F)
cbg2 <- filter(cbg2, true == T)
unique(cbg2$Species) #165 species


## Talking to Cheryl
cbg <- read.csv("Projects/Seed-Traits/Datasets/RSABG-Dryland-annuals_2021.csv") 
cbg.edit <- read.csv("Projects/Seed-Traits/Datasets/Orders/CBG-Cheryl_edit.csv")

cbg <- filter(cbg, !Species %in% cbg.edit$Species)

cbg.stat <- filter(cbg.stat, CRPR != "1B.1", CRPR != "1B.2", CRPR != "4", CRPR != "1B.3", CRPR != "1B", CRPR != "2B.2" , CRPR != "2" , CRPR != "4.3", CRPR != "2B.1", CRPR != "4.2", CRPR != "2B.3", CRPR != "2.3", CRPR != "3", CRPR != "3.2", CRPR != "2.2", CRPR != "2B")

cbg <- filter(cbg, Species %in% cbg.stat$NAME, date.store > 1995)
cbg$seed.num.1 <- cbg$seed.num*0.1
cbg$seed.15g <- cbg$seed.per.g*0.15
write.csv(cbg, "Projects/Seed-Traits/Datasets/Orders/CBG-add.csv", row.names = F)

grin.new <- read.csv("Projects/Seed-Traits/Datasets/Orders/grin3.csv")
grin <- filter(grin, Species %in% grin.new$Species)
write.csv(grin, "Projects/Seed-Traits/Datasets/Orders/grin-order.csv", row.names = F)

####
grin <- read.csv("Projects/Seed-Traits/Datasets/grin-order-seedmass.csv")
df <- read.csv("Projects/Seed-Traits/Datasets/Orders/Seed_Obtained.csv")
df <- merge(df, grin, by = "Species")