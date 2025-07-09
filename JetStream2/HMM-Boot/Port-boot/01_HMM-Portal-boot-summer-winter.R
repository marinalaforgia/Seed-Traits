rm(list=ls())

source("~/Documents/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(doParallel)
library(portalr)
library(foreach)

plant_data <- load_plant_data(
  path = "/home/exouser/Documents/Portal-boot/",
  download_if_missing = T,
  quiet = FALSE
)

plant.data <- clean_plant_data(
  plant_data,
  type = 'annual',
  unknowns = FALSE,
  correct_sp = TRUE
)


plant.data <- filter(plant.data, community == "Summer and Winter Annual")

# treatments change over time (See https://github.com/weecology/PortalData/blob/main/SiteandMethods/Portal_plot_treatments.csv) so although the project dates back to 1977, running models on controls must be split into different segments:

controls <- c(2,11,14,22) #controls from 1988 - 2009; 2010 not censused

plant.data <- filter(plant.data, year >= 1988 & year < 2010, plot %in% controls)


#### Summer AND winter ####
#### .___Prep Data ####
portal <- plant.data

portal$PA <- 1

# for summer AND winter species (s1 = season 1 = winter, s2 = season 2 = summer; winter census always occurs in march/april and summer census always occurs in august/september)
portal$year.season <- ifelse(portal$season == "summer", paste(portal$year, "s2", sep = "."), paste(portal$year, "s1", sep = "."))
portal <- portal[,c(19,3:5,20)]
portal <- portal %>%
  dplyr::group_by(year.season, plot, quadrat, species) %>%
  dplyr::summarize(PA = mean(PA))

portal <- pivot_wider(portal, names_from = "year.season", values_from = "PA")

portal[is.na(portal)] <- 0

# according to the metadata these were all done https://github.com/weecology/PortalData/blob/main/Plants/Portal_plant_censuses.csv
portal$'1996.s1' <- 0   
portal$'2000.s1' <- 0  
portal$'2009.s2' <- 0  
portal$'2009.s1' <- 0  
portal$'2003.s2' <- 0  
portal$'2006.s1' <- 0  


portal <- portal[ , order(names(portal))]
portal  <- portal[,c(45:47, 1:44)]
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

port.df <- expand.grid(Species_Name = names(species.list), boot = 1:100, p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA, n.plots = NA, min.years = NA, max.years = NA, mean.years = NA) # empty df to fill with rates

#Setup backend to use many processors
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-2) 
registerDoParallel(cluster)

boot.df <- foreach(j = names(species.list), .combine = rbind) %dopar% { # for each species
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
    
    port.df[port.df$Species_Name == j & port.df$boot == k,]$p0 <- EMresult$param$p0
    port.df[port.df$Species_Name == j & port.df$boot == k,]$g <- EMresult$param$g
    port.df[port.df$Species_Name == j & port.df$boot == k,]$c <- EMresult$param$c
    port.df[port.df$Species_Name == j & port.df$boot == k,]$s <- EMresult$param$s
    port.df[port.df$Species_Name == j & port.df$boot == k,]$r <- EMresult$param$r
    port.df[port.df$Species_Name == j & port.df$boot == k,]$iter <- which.max(EMresult$lllist)
    
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

colnames(boot.df) <- c("Species_Name", "p0.boot", "g.boot", "c.boot", "s.boot", "r.boot", "iter.boot","p0.boot.sd", "g.boot.sd", "c.boot.sd", "s.boot.sd", "r.boot.sd", "iter.boot.sd", "perc.con")


saveRDS(boot.df, "port-species-boot-summer-winter.RDS")
