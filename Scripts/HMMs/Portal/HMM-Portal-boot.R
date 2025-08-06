rm(list=ls())
source("Scripts/HMMs/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)

plant_data <- load_plant_data(
  path = "Data/Long-Term-Datasets/Portal/",
  download_if_missing = F,
  quiet = FALSE
)

plant.data <- clean_plant_data(
  plant_data,
  type = 'annual',
  unknowns = FALSE,
  correct_sp = TRUE
)


plant.data <- filter(plant.data, season == "winter", community == "Winter Annual") # only winter census and only winter SPECIES (no species that occur in both seasons as this will affect P/A but since we aren't looking at summer and they could be there in summer these species' estiamtes may be incorrect

# treatments change over time (See https://github.com/weecology/PortalData/blob/main/SiteandMethods/Portal_plot_treatments.csv) so although the project dates back to 1977, running models on controls must be split into different segments:

controls <- c(2,11,14,22) #controls from 1988 - 2009; 2010 not censused

plant.data <- filter(plant.data, year >= 1988 & year < 2010, plot %in% controls)


#### Winter ####
#### .___Prep Data ####
portal <- plant.data

portal$PA <- 1
portal <- portal[,c(1,3:5,19)]
portal <- portal %>%
  dplyr::group_by(year, plot, quadrat, species) %>%
  dplyr::summarize(PA = mean(PA))

portal <- pivot_wider(portal, names_from = "year", values_from = "PA")

portal[is.na(portal)] <- 0

# # according to the metadata 1996, 2000, 2009 were all done... https://github.com/weecology/PortalData/blob/main/Plants/Portal_plant_censuses.csv
portal$'1996' <- 0   #winter
portal$'2000' <- 0  # winter
portal$'2009' <- 0  # winter


portal <- portal[ , order(names(portal))]
portal <- portal[,c(23:25, 1:22)]
colnames(portal)[3] <- "SpeciesCode"

#### .___Create Species List ####

species.list <- list()

for(i in unique(portal$SpeciesCode)) {
  tmp <- filter(portal, SpeciesCode == i)
  tmp <- tmp[,-c(1:3)]
  tmp <- tmp[ , sort(colnames(tmp))]
  species.list[[i]] <- tmp 
}


#### Run HMM ####

port.df <- expand.grid(Species_Name = names(species.list), sim = 1:100, p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA, n.plots = NA, min.years = NA, max.years = NA, mean.years = NA) # empty df to fill with rates

#Setup backend to use many processors
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-1) 
registerDoParallel(cluster)

sim.df <- foreach(j = names(species.list)[1], .combine = rbind) %dopar% { # for each species
  for(k in 1:3) {
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
  
  port.df[port.df$Species_Name == j & port.df$sim == k,]$p0 <- EMresult$param$p0
  port.df[port.df$Species_Name == j & port.df$sim == k,]$g <- EMresult$param$g
  port.df[port.df$Species_Name == j & port.df$sim == k,]$c <- EMresult$param$c
  port.df[port.df$Species_Name == j & port.df$sim == k,]$s <- EMresult$param$s
  port.df[port.df$Species_Name == j & port.df$sim == k,]$r <- EMresult$param$r
  port.df[port.df$Species_Name == j & port.df$sim == k,]$iter <- which.max(EMresult$lllist)
  
  }

      cbind(j, 
            mean(port.df[port.df$Species_Name == j,]$p0), 
            mean(port.df[port.df$Species_Name == j,]$g), 
            mean(port.df[port.df$Species_Name == j,]$c), 
            mean(port.df[port.df$Species_Name == j,]$s),
            mean(port.df[port.df$Species_Name == j,]$r),
            mean(port.df[port.df$Species_Name == j,]$iter),
            sd(port.df[port.df$Species_Name == j,]$p0), 
            sd(port.df[port.df$Species_Name == j,]$g), 
            sd(port.df[port.df$Species_Name == j,]$c), 
            sd(port.df[port.df$Species_Name == j,]$s),
            sd(port.df[port.df$Species_Name == j,]$r),
            sd(port.df[port.df$Species_Name == j,]$iter),
            nrow(port.df[which(port.df$Species_Name == j & port.df$iter < 150),])/100
      )
}

#Stop cluster
stopCluster(cluster)

colnames(sim.df) <- c("Species_Name", "p0.boot", "g.boot", "c.boot", "s.boot", "r.boot", "iter.boot","p0.boot.sd", "g.boot.sd", "c.boot.sd", "s.boot.sd", "r.boot.sd", "iter.boot.sd", "perc.con")

sim.df <- as.data.frame(sim.df)
saveRDS(sim.df, "port-species-boot.RDS")
