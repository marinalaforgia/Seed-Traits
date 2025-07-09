# Annuals across different lists
library(plyr)
library(tidyverse)

RSABG.old <- read.csv("Seed-Databases/RSABG-Online-Seed-List-Jan-2020.csv")
Calflora <- read.csv("Seed-Databases/Calflora_annuals-w-families.csv")
Calflora.DA <- read.csv("Seed-Databases/Cal-Flora_dryland-annuals.csv")
#ProjBase <- read.csv("Seed-Databases/Project-Baseline_Species.csv")
SBBG <- read.csv("Seed-Databases/20191101_SBBG_Seed-Bank.csv")
DBG <- read.csv("Seed-Databases/Collection_Data_Results_DBG.csv")
AZ.annuals <- read.csv("Seed-Databases/Arizonia_Annuals.csv") 
chosen <- read.csv("Seed-Databases/20210112_Possible-Species.csv")

#### RSABG ####
# Filter to just Annuals with more than 250 seeds in the collection
#RSABG <- filter(RSABG, NAME %in% Calflora$Species, PT != "G", STATUS == " ", YEAR > 2000) 
RSABG <- filter(RSABG, name %in% Calflora$Species) 
write.csv(RSABG, "RSABG-annuals-tmp-20210121.csv", row.names = F)
RSABG <- merge(RSABG, Calflora, by.x = "NAME", by.y = "Species", all.y = F)
# RSABG.new <- as.data.frame(matrix(NA, 361, 21))
# colnames(RSABG.new) <- colnames(RSABG)

# test <- ddply(RSABG, .(scientificName, Family), summarize, organismQuantitynumberOfSeeds[which(organismQuantitynumberOfSeeds == max(organismQuantitynumberOfSeeds))])
# colnames(test)[3] <- "organismQuantitynumberOfSeeds"
# 
# test2 <- merge(test, RSABG, by = c("scientificName","Family", "organismQuantitynumberOfSeeds"), all.y = F)

#test3 <- sample_n(test2, size = 150, replace = F)

# for(i in 1:length(RSABG.sp)){
#   RSABG.new[i,] <- RSABG[which.max(as.character(RSABG$scientificName) == as.character(RSABG.sp[i]) & RSABG$organismQuantitynumberOfSeeds),]
# }


write.csv(test2, "RSABG-Annuals-no-duplicates.csv", row.names = F)
write.csv(RSABG, "RSABG-Annuals.csv", row.names = F)

# RSABG$eventDate <- as.Date(RSABG$eventDate, "%d %B %Y")
# write.csv(RSABG, "RSABG-Annuals.csv", row.names = F)
RSABG.sp <- unique(as.character(RSABG$scientificName))

# ProjBase <- filter(ProjBase, longevity == "Annual") 
# PB.sp <- unique(ProjBase$name) # 12 Species

#### SBBG ####
SBBG$Scientific.name <- gsub("\\s*\\([^\\)]+\\)", "", as.character(SBBG$Scientific.name))
not <- c("USFS", "US Navy", "TNC and Energy plant", "U.S. Navy", "National Park Service")
SBBG <- filter(SBBG, Scientific.name %in% Calflora$Species, Number.of.Propagules > 500, !Landowner %in% not, Ownership.of.Seeds == "SBBG", Fed.Status == "", State.Status == "") 
SBBG.sp <- unique(as.character(SBBG$Scientific.name)) # 23 Species

write.csv(SBBG, "SBBG-annuals.csv", row.names = F)
#### DBG ####
DBG <- filter(DBG, State.Province == "Arizona")

DBG$Name <- trimws(DBG$Name, which = "right")

# Filter to just Annuals 
DBG <- filter(DBG, Name %in% AZ.annuals$Scientific.Name) 
test <- unique(DBG$Name)

write.csv(DBG, "DBG-annuals.csv", row.names = F)

#### Kew ####
kew <- read.csv("Seed-Databases/Kew-Millenium-AZ.csv")
kew$ScientificName <- paste(kew$Genus, kew$Species, kew$Subsp., kew$Subsp2, sep = " ")
kew$ScientificName <- trimws(kew$ScientificName, which = "right")
kew$ScientificName <- filter(kew, ScientificName %in% AZ.annuals$Scientific.Name) 
kew <- filter(kew, ScientificName != "Boerhavia erecta", ScientificName != "Penstemon pachyphyllus", ScientificName != "Oenothera pubescens")

write.csv(kew, "kew-annuals-AZ.csv", row.names = F)

SD.species <- read.csv("Sonoran-Desert_Species-List_2016.csv")
SD.species$Name <- paste(SD.species$Genus, SD.species$Species, sep = " ")

DBG.herb <- read.csv("Seed-Databases/DBG_Herbarium_occurrences.csv")
DBG.herb <- filter(DBG.herb, scientificName %in% AZ.annuals$Scientific.Name) 

#DBG.herb1 <- filter(DBG.herb, scientificName %in% SD.species$Name | scientificName %in% kew$ScientificName, year > 2000) # 98 species

#DBG.herb3 <- filter(DBG.herb, scientificName %in% SD.species$Name)
DBG.herb <- filter(DBG.herb, scientificName %in% kew$ScientificName, year > 2000) 

#DBG.herb <- filter(DBG.herb, scientificName %in% SD.species$Name) 
#unique(DBG.herb1$scientificName)
#DBG.herb2 <- DBG.herb1 %>% group_by(scientificName) %>% sample_n(1)
write.csv(DBG.herb, "DBG-herbaria-kew.csv", row.names = F)


#### GRIN ####
grin <- read.csv("Post-Doc/NSF-PRFB/Seed-Databases/USDA-GRIN.csv")
grin$Received <- as.numeric(as.character(grin$Received))
grin <- filter(grin, Species %in% AZ.annuals$Scientific.Name, Improvement.Status == "WILD", Received > 1990)
write.csv(grin, "USDA-GRIN-annuals.csv", row.names = F)
grin.sp <- unique(grin$Species)
DBG[DBG$Name %in% grin$Species,]$Name
# 64 species in GRIN plus 22 DBG + 

#### FAMILIES ####
families <- c(as.character(Calflora$Family), as.character(AZ.annuals$Family))
families <- unique(families) # 80 families

grin <- merge(grin, AZ.annuals[,2:3], by.x = "Species", by.y = "Scientific.Name", all.x = T)
grin.fam <- unique(grin$Family)

DBG.fam <- unique(DBG$Family)

RSABG.fam <- unique(RSABG$Family)

SBBG <- merge(SBBG, Calflora, by.x = "Scientific.name", by.y = "Species", all.y = F)

SBBG.fam <- unique(SBBG$Family)

families.sam <- c(as.character(SBBG.fam), as.character(RSABG.fam), as.character(grin.fam), as.character(DBG.fam))
families.sam <- unique(families.sam) # 71 families

# Chosen species
chosen <- merge(chosen, Calflora, by = "Species", all.x = T)
chosen <- merge(chosen, AZ.annuals[,c(2,3)], by.x = c("Species", "Family"), by.y = c("Scientific.Name", "Family"), all.x = T)

RSABG.sp <- unique(RSABG[,c("NAME", "Family")])
RSABG.sp$Site <- "RSABG"
RSABG.sp$FunGroup <- NA
RSABG.sp$Nat.Ex <- NA
RSABG.sp$Type <- NA
colnames(RSABG.sp)[1] <- "Species"
chosen <- rbind(chosen, RSABG.sp)
chosen[which(duplicated(chosen$Species) == T),]

grin.sp <- unique(grin[,c("Species", "Family")])
grin.sp$Site <- "GRIN"
grin.sp$FunGroup <- NA
grin.sp$Nat.Ex <- NA
grin.sp$Type <- NA
chosen <- rbind(chosen, grin.sp)

DBG.sp <- unique(DBG[,c("Name", "Family")])
colnames(DBG.sp)[1] <- "Species"
DBG.sp$Site <- "DBG"
DBG.sp$FunGroup <- NA
DBG.sp$Nat.Ex <- NA
DBG.sp$Type <- NA
chosen <- rbind(chosen, DBG.sp)

#----------
RSABG <- read.csv("RSABG-annuals-tmp-20210121.csv")
RSABG.old <- filter(RSABG.old, STATUS != " ", STATUS != "R") 
RSABG <- filter(RSABG, !RSABG$name %in% RSABG.old$NAME)
RSABG$seed_.5g <- RSABG$seed.per.g*0.5
RSABG <- filter(RSABG, 0.1*seed.num > seed_.5g, date.store > 2000)
unique(RSABG$name) #112 species, plus ~20 from other stuff = 
