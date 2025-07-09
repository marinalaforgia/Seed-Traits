library(plyr)
library(dplyr)
library(tidyverse)
library(openxlsx)

rm(list=ls())

# setwd("/Users/Marina/Documents/USDA-PostDoc/Projects/Seed-Traits/Data/Shape_Size/Videometer/")
# filenames <- list.files("raw", pattern="*.xlsx", full.names=TRUE)
# 
# full <- data.frame()
# 
# for(i in 1:length(filenames)) {
#   tmp <- read_excel(filenames[i])
#   full <- rbind(full, tmp)
# }
# 
# full.test <- full %>% separate(Org.Filename, c("Species.Code", "ID", "chem.morph", "2d.3d", "seg.var", "extra1", "extra2"))
# 
# 
# write.csv(full.test, "20220221_Full.csv", row.names = F)

# Insert 4 hours of dumb hand cleanup because file name metadata were not in order and pieces were missing

full <- read.csv("Data/Shape_Size/Archive/20220221_Videometer_raw.csv")

vid.ch <- read.xlsx("Data/Shape_Size/Videometer-Process/Videometer-checklist.xlsx")

# chem.morph.pic = unit taken a pic of and processed comparable to videometer check list
# chem.morph.fun = function unit comparable to all seeds accessions and mass
full.test <- merge(full, vid.ch, by.x = c("Species.Code","ID", "chem.morph.pic"), by.y = c("code", "ID.short", "Unit"), all = T) 

full.test1 <- anti_join(full, vid.ch, c("Species.Code" = "code","ID" = "ID.short","chem.morph.pic" = "Unit")) # none

full.test2 <- anti_join(vid.ch, full, c("code" = "Species.Code","ID.short" = "ID", "Unit" = "chem.morph.pic")) # ones we havent done yet



full.2d <- filter(full, X2d.3d == "2D")
full.test.2d <- merge(full.2d, vid.ch, by.x = c("Species.Code","ID","chem.morph.pic"), by.y = c("code", "ID.short", "Unit"), all=T) ## Everything is there!

full.3d <- filter(full, X2d.3d == "3D")
full.test.3d <- merge(full.3d, vid.ch[vid.ch$X3d == "y",], by.x = c("Species.Code","ID","chem.morph.pic"), by.y = c("code", "ID.short", "Unit"), all=T) ## Everything is there!

full.test <- full.test[,c(33,1,35,20:22,5:10)]
colnames(full.test)[3] <- "ID"
colnames(full.test)[6] <- "app"
write.csv(full.test, "Data/Shape_Size/20220221_Videometer.csv", row.names = F)
