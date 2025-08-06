#### Carrizo HMM ####
rm(list=ls())
source("Scripts/HMMs/HMM_Functions.R")
# Swain is more species rich and dif soil type from center

#### Load Libraries ####
library(ggplot2)
library(tidyverse)

#### Prep Data ####
cp <- read.csv("Data/Long-Term-Datasets/Carrizo-Plain/Carrizo_2007_2021_coverdata_long_format_with_cover_percentage_final.csv")

# LASGRA accidentally input as LASGLA in 2017
cp[cp$SpeciesCode == "LASGLA",]$SpeciesCode <- "LASGRA"

sort(unique(cp$SpeciesCode))

rm <- c("BARE", "LITTER","FRESHDIRT","BPHOLE","HOLE","ASTspp","COWPIE","VULPIA","ANT","ALLIUM","BURROW","GOPHER","unknown10", "UNK 10", "UNK 11", "MOSS", "DICCAP")

### Filter out sites with less than 15 years of data (these are sites with gaps and gets rid of sites that moved for the precip exp)
cp.sum <- cp %>%
  dplyr::group_by(New.Plot.ID) %>%
  dplyr::summarize(n.years = length(unique(year))) %>%
  filter(n.years == 15)# Some plots were moved in 2015 to better align with the precipitation experiment. as long as I keep the entire time series (15 years), I can include controls as those that have data for the entire time series did not move

precip <- c("CONTROL", "NONE")

#removing the first year of data because 75% of plots are 0s for everything, this doesnt seem to affect center as much as it seems to affect Swain (but swain like 90% are zeros and with such a smaller dataset it's worse)
cp <- filter(cp, New.Plot.ID %in% cp.sum$New.Plot.ID, !(SpeciesCode %in% rm), cattle == "Ungrazed", GKR == "GKR", Precip_Treatment %in% precip, year > 2007)

#cp <- filter(cp, New.Plot.ID %in% cp.sum$New.Plot.ID, !(SpeciesCode %in% rm), GKR == "Excl", Precip_Treatment %in% precip, year > 2007)

#for Hannah
#write.csv(cp, "20230907_CP-data-for-Hanna.csv", row.names = F)

length(unique(cp$New.Plot.ID))

## Relative abundance for cwm traits ###
cp.RA <- cp %>%
  dplyr::mutate(sum.abundance = sum(Count2)) %>%
  dplyr::group_by(SpeciesCode, sum.abundance) %>%
  dplyr::summarize(sp.RA = sum(Count2)) %>%
  dplyr::group_by(SpeciesCode) %>%
  dplyr::summarize(RA = sp.RA/sum.abundance)

# dups <- cp %>%
#     dplyr::group_by(New.Plot.ID, SpeciesCode, year) %>%
#     dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
#     dplyr::filter(n > 1L) 

cp$PA <- 1
cp <- pivot_wider(cp[,c(2,3,11,22)], names_from = "year", values_from = "PA")
cp[is.na(cp)] <- 0
cp[cp$SpeciesCode == "PHLGRA",]$SpeciesCode <- "MICGRA" # update species name

#### Metadata ####
cp$sumPA <- rowSums(cp[,3:16])

meta.cp <- cp %>%
  dplyr::group_by(SpeciesCode) %>%
  dplyr::summarize(n.plots = length(New.Plot.ID),
            mean.years = mean(sumPA),
            min.years = min(sumPA),
            max.years = max(sumPA))

cp <- cp[,-17]

#### Create Species List ####
length(unique(cp$New.Plot.ID))

species.list <- list()

cp <- filter(cp, SpeciesCode != "EROCIC") # code breaks with Erodium cicuatrium because it's in too many plots in too many years

for(i in unique(cp$SpeciesCode)) {
  tmp <- filter(cp, SpeciesCode == i)
  tmp <- tmp[,-c(1:2)]
  tmp <- tmp[ , sort(colnames(tmp))]
  species.list[[i]] <- tmp 
}


#### Run HMM ####

CP.df <- expand.grid(Species_Name = names(species.list), p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA, n.plots = NA, min.years = NA, max.years = NA, mean.years = NA) # empty df to fill with rates

for(j in names(species.list)){ # for each species
  X = as.matrix(species.list[[j]]) # use their time series PA data
  n = 5 # why n = 5?
  p0Results = rep(0,n)
  gResults = rep(0,n)
  cResults = rep(0,n)
  sResults = rep(0,n)
  rResults = rep(0,n)
  llResults = rep(0,n)
  
print(j)
print(j)  
print(j)  
print(j)  
print(j)  

  for (i in 1:n) {
    for (k in 1:5) print(i) 
    EMresult = EMestimation(X, r = 1) # update the model params, stop when new log likelihood is lower than the old log likelihood plus some precision 
    # fill in results with new params from the best log likelihood model
    p0Results[i] = EMresult$param$p0 # initial seed bank prob
    gResults[i] = EMresult$param$g # germ
    cResults[i] = EMresult$param$c #col
    sResults[i] = EMresult$param$s #surv
    rResults[i] = EMresult$param$r #prod
    llResults[i] = EMresult$ll # log-likelihood

  }

  CP.df[CP.df$Species_Name == j,]$p0 <- EMresult$param$p0
  CP.df[CP.df$Species_Name == j,]$g <- EMresult$param$g
  CP.df[CP.df$Species_Name == j,]$c <- EMresult$param$c
  CP.df[CP.df$Species_Name == j,]$s <- EMresult$param$s
  CP.df[CP.df$Species_Name == j,]$r <- EMresult$param$r
  CP.df[CP.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
  CP.df[CP.df$Species_Name == j,]$n.plots <- nrow(species.list[[j]])
  CP.df[CP.df$Species_Name == j,]$mean.years <- mean(rowSums(species.list[[j]]))
  CP.df[CP.df$Species_Name == j,]$min.years <- min(rowSums(species.list[[j]]))
  CP.df[CP.df$Species_Name == j,]$max.years <- max(rowSums(species.list[[j]]))
  
}

#### Explore Output ####
trait <- read.csv("/Users/Marina/Documents/USDA-PostDoc/Projects/Seed-Traits/Data/Long-Term-Datasets/Carrizo-Plain/Carrizo_Traits.csv")
trait$FunGroup <- paste(trait$Status, trait$Habit)
colnames(CP.df)[1] <- "Code"
CP.df <- merge(CP.df, unique(trait[,c(1,2,10,12)]), by.x = "Code", by.y = "SpeciesCode", all.x = T)
CP.df <- filter(CP.df, Annual.Perennial == "Annual")
CP.df <- CP.df[,c(12,1,14,2:11)]

meta.cp <- merge(meta.cp, unique(trait[,c(1,2,10,12)]), by = "SpeciesCode", all.x = T)
meta.cp <- filter(meta.cp, Annual.Perennial == "Annual")
meta.cp <- meta.cp[,c(6,1,8,2,4,5,3)]
write.csv(meta.cp, "Data/Long-Term-Datasets/Carrizo-Plain/HMM-meta-CP.csv", row.names = F)


# ### Save RA data ###
# saveRDS(cp.RA, "Data/Long-Term-Datasets/Carrizo-Plain/Carrizo_RA.RDS")

ggplot(CP.df[CP.df$iter < 150 & CP.df$n.plots >= 15,], aes(x = s, y = c)) +
  #geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  geom_text(aes(label = Code)) +
  theme_classic()


saveRDS(CP.df, "Scripts/HMMs/Carrizo-Plain/CP_HMM_none-control_ungrazed_GKR_2008.RDS")

