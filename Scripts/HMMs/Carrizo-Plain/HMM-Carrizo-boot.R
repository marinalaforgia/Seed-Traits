#### Carrizo Plain HMM Boot ####

rm(list=ls())
source("HMM_Functions.R")
# Swain is more species rich and dif soil type from center

#### Load Libraries ####
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(doParallel)
library(portalr)
library(foreach)

#### Prep Data ####
cp <- read.csv("Carrizo_2007_2021_coverdata_long_format_with_cover_percentage_final.csv")

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


length(unique(cp$New.Plot.ID))

cp$PA <- 1
cp <- pivot_wider(cp[,c(2,3,11,22)], names_from = "year", values_from = "PA")
cp[is.na(cp)] <- 0
cp[cp$SpeciesCode == "PHLGRA",]$SpeciesCode <- "MICGRA" # update species name

#### Create Species List ####
length(unique(cp$New.Plot.ID))


species.list <- list()

cp <- filter(cp, SpeciesCode != "EROCIC") # code breaks with Erodium cicuatrium because it's in too many plots in too many years

#### Create Species List ####


species.list <- list()

for(i in unique(cp$SpeciesCode)) {
  tmp <- filter(cp, SpeciesCode == i)
  tmp <- tmp[,-c(1:2)]
  tmp <- tmp[ , sort(colnames(tmp))]
  species.list[[i]] <- tmp 
}

#### Run HMM ####

CP.df <- expand.grid(Species_Name = names(species.list), sim = 1:100, p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA, n.plots = NA, min.years = NA, max.years = NA, mean.years = NA) # empty df to fill with rates

#Setup backend to use many processors
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-4) 
registerDoParallel(cluster)

sim.df <- foreach(j = names(species.list), .combine = rbind) %dopar% { # for each species
  for(k in 1:100) {
  X = as.matrix(species.list[[j]]) # use their time series PA data
  n = 5 
  p0Results = rep(0,n)
  gResults = rep(0,n)
  cResults = rep(0,n)
  sResults = rep(0,n)
  rResults = rep(0,n)
  llResults = rep(0,n)
 
  print(j)  
  
  for (i in 1:n) {
    for (a in 1:5) print(i) # I don't understand what's going on here, why print it 5 times?
    EMresult = EMestimation(X, r = 1) # update the model params, stop when new log likelihood is lower than the old log likelihood plus some precision 
   
    # fill in results with new params from the best log likelihood model
    p0Results[i] = EMresult$param$p0 # initial seed bank prob
    gResults[i] = EMresult$param$g # germ
    cResults[i] = EMresult$param$c #col
    sResults[i] = EMresult$param$s #surv
    rResults[i] = EMresult$param$r #prod
    llResults[i] = EMresult$ll # log-likelihood

  }
  
  CP.df[CP.df$Species_Name == j & CP.df$sim == k,]$p0 <- EMresult$param$p0
  CP.df[CP.df$Species_Name == j & CP.df$sim == k,]$g <- EMresult$param$g
  CP.df[CP.df$Species_Name == j & CP.df$sim == k,]$c <- EMresult$param$c
  CP.df[CP.df$Species_Name == j & CP.df$sim == k,]$s <- EMresult$param$s
  CP.df[CP.df$Species_Name == j & CP.df$sim == k,]$r <- EMresult$param$r
  CP.df[CP.df$Species_Name == j & CP.df$sim == k,]$iter <- which.max(EMresult$lllist)
  
  }

      cbind(j, 
            mean(CP.df[CP.df$Species_Name == j,]$p0), 
            mean(CP.df[CP.df$Species_Name == j,]$g), 
            mean(CP.df[CP.df$Species_Name == j,]$c), 
            mean(CP.df[CP.df$Species_Name == j,]$s),
            mean(CP.df[CP.df$Species_Name == j,]$r),
            mean(CP.df[CP.df$Species_Name == j,]$iter),
            sd(CP.df[CP.df$Species_Name == j,]$p0), 
            sd(CP.df[CP.df$Species_Name == j,]$g), 
            sd(CP.df[CP.df$Species_Name == j,]$c), 
            sd(CP.df[CP.df$Species_Name == j,]$s),
            sd(CP.df[CP.df$Species_Name == j,]$r),
            sd(CP.df[CP.df$Species_Name == j,]$iter,
            nrow(CP.df[CP.df$Species_Name == j & CP.df$iter < 150,])/100)
      )
}

#Stop cluster
stopCluster(cluster)

colnames(sim.df) <- c("Species_Name", "p0.boot", "g.boot", "c.boot", "s.boot", "r.boot", "iter.boot","p0.boot.sd", "g.boot.sd", "c.boot.sd", "s.boot.sd", "r.boot.sd", "iter.boot.sd")


saveRDS(sim.df, "CP-species-boot.RDS")
