#### Herbarium records ####

# script that takes herbarium record occurences for all of my seed trait species to calculate
# mean AI
# min AI
# max AI
# quartile AI
# standard deviation AI

rm(list=ls())

#### Load Libraries ####
library(raster)
library(rgdal)
library(tidyverse)
library(terra)
library(leaflet)
library(maps)
library(sp)
library(sf)
library(Taxonstand)
library(CoordinateCleaner)

####  Load Data ####
# herbarium records for all species that i have seed trait data for
traits <- read.csv("Data/20230605_Seed-Traits_clean.csv")
occ.cch2 <- read.csv("Data/Herbarium-Occurences/CCH2/occurrences.csv") 
occ.sein <- read.csv("Data/Herbarium-Occurences/SEIN/occurrences.csv") 

occ <- rbind(occ.cch2, occ.sein)

occ <- occ[!is.na(occ$decimalLatitude),]
occ <- occ[!is.na(occ$decimalLongitude),]

#### Fix herbarium names ####
# convert species names to names in my species list
length(unique(traits$Species)) #278 species 
length(unique(occ$scientificName)) #2073 species output for 278 species input

occ.names <- unique(occ$scientificName)

tpl.names <- TPL(occ.names)

tpl.names$tpl.name <- trimws(paste(tpl.names$New.Genus, tpl.names$New.Species, tpl.names$New.Infraspecific.rank, tpl.names$New.Infraspecific, sep = " "))

colnames(tpl.names)[1] <- "herb.name"
tpl.names <- tpl.names[,c(1,26,2:25)]

#by.hand <- unique(tpl.names$tpl.name) # 632 unique names

#write.csv(by.hand, "Spatial-Data/herbarium-name-fix.csv", row.names = F)

# had to go through by hand and correct these names/find synonyms
by.hand <- read.csv("Data/Herbarium-Occurences/herbarium-name-fix.csv")

tpl.names <- merge(by.hand, tpl.names, by = "tpl.name", all.x = F, all.y = T)

occ <- merge(tpl.names[,c(2,3)], occ, by.x = "herb.name", by.y = "scientificName") # 278 unique species!! YAY

#### Clean geographic outliers ####
#https://search.r-project.org/CRAN/refmans/CoordinateCleaner/html/cc_outl.html
# remove duplicates
occ <- cc_dupl(
  occ,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  species = "Marina.species") 

# one specimen has a latitude of 334 and another a longitude over 180
occ <- filter(occ, decimalLatitude <= 90)
occ <- filter(occ, decimalLongitude <= 180)

# remove specimen outliers using quantile method
occ <- cc_outl(occ,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  species = "Marina.species") #removed  4638 outliers 


#### Extract aridity data ####
# load in entire map of aridity data
r <- raster("Spatial-Data/ai_et0/ai_et0.tif")

# cut down to US and mexico to

countries <-rgdal::readOGR("Spatial-Data/shapefiles/world-administrative-boundaries/world-administrative-boundaries.shp")

countries <- countries[countries$name == "United States of America"|
                       countries$name == "Mexico"
                        ,]

r <- crop(r, countries)

# spatial resolution of AI data is about 30 arc seconds (roughly 1 km), so creating a buffer around each point at the same level 
occ_buffer <- st_buffer(st_as_sf(occ, coords = c("decimalLongitude", "decimalLatitude"), crs = "WGS84"), 30) #create buffer of 30 arc degrees around each point

AI <- terra::extract(r, occ_buffer, mean)

occ <- cbind(occ, AI)

saveRDS(occ, "Data/Herbarium-Occurences/Herbarium_occ_AI.RDS")

occ <- readRDS("Data/Herbarium-Occurences/Herbarium_occ_AI.RDS")
#### Calculate AI variables ####
# mean AI
# min AI
# max AI
# quartile AI
# standard deviation AI

occ$AI <- occ$AI*0.0001

occ <- filter(occ, decimalLatitude > 0.1) # gets rid of erroneous coordinate for a CA specimen recorded as 0.000034
occ <- filter(occ, decimalLatitude <= 49) # gets rid of things north of canada border
occ <- filter(occ, stateProvince != "Hawaii") 
occ <- filter(occ, decimalLongitude > -127.6522) # two coordinates in the middle of the ocean

occ.CESO <- filter(occ, Marina.species == "Centaurea solstitialis")

leaflet(occ.CESO) %>%
  addTiles() %>%
  addRasterImage(r) %>%
  addMarkers(lng = occ.CESO$decimalLongitude, 
             lat = occ.CESO$decimalLatitude, 
             label = occ.CESO$ID)

occ.sum <- occ %>%
  group_by(Marina.species) %>%
  summarise(AI.IQR = IQR(AI, na.rm = T),
            AI.25 = quantile(AI, na.rm = T)[2],
            AI.mean = mean(AI, na.rm = T),
            AI.75 = quantile(AI, na.rm = T)[4]
            )

saveRDS(occ.sum, "Data/Herbarium-Occurences/Herbarium_occ_AI_sum.RDS")

#### Explore ####
occ <- readRDS("Data/Herbarium-Occurences/Herbarium_occ_AI.RDS")

#### Calculate Annual AI at each site ####
r <- raster("Spatial-Data/Global-AI_ET0_v3_annual/ai_v3_yr.tif")
