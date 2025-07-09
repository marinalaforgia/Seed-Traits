### Add AI to Accessions 05/30/2023 ###
rm(list=ls())

library(tidyverse)
library(terra)
library(leaflet)
library(maps)
library(raster)
library(sf)

acc <- read.csv("Data/20230530_Seeds_All-Accessions.csv") # full final accession csv
acc2 <- filter(acc, !is.na(lat))

r <- raster("Spatial-Data/AI-Raster.tif")

leaflet(sp_buffer) %>%
  addTiles() %>%
  addRasterImage(r) %>%
  addPolygons() %>% 
  addMarkers(lng = sp_buffer$long, 
             lat = sp_buffer$lat, 
             label = sp_buffer$ID)

sp_buffer <- st_buffer(st_as_sf(acc2, coords = c("long", "lat"), crs = "WGS84"), 100) #create buffer of 3 arc seconds around point to calculate AI, this likely will do nothing for those inside a 30 arc second pixel (resolution of AI data), but should help for those points right on the border between two pixels

AI <- terra::extract(r, sp_buffer, mean)

acc2 <- cbind(acc2, AI)

acc2$AI <- acc2$AI*0.0001

acc <- merge(acc, acc2[,c(1,30)], by = "ID", all = T)

write.csv(acc, "Data/20230530_Seeds_All-Accessions.csv", row.names = F)
