#### Sonoran Desert HMM Boot ####

rm(list=ls())
source("HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(doParallel)
library(portalr)
library(foreach)

#### Prep Data ####
sd <- read.csv("SD-boot/CensusData.csv")

iffy.sites <- c("3-open-west", "3-open-east", "0/3-shrub-north", "30-open-north", "30-open-south", "30-shrub-north", "30-shrub-south", "7-open-east", "7-open-west", "8-open-east", "8-open-west")

#remove weird sites
sd <- filter(sd, !(plot.habitat.replicate %in% iffy.sites))
length(unique(sd$plot.habitat.replicate))
plot.size <- read.csv("SD-boot/plot-sizes.csv")


remove <- c("CRsp", # removing unidentified species and known perennials
            "DA-SP",
            "eriog",
            "Eriog",
            "ERIOG",
            "LO-AS",
            "LOsp",
            "PHsp",
            "PLsp",
            "plsp",
            "ersp",
            "crsp",
            "losp"
            )

trait <- read.csv("SD-boot/Species_list_2016_2.csv")
trait$Species_Name <- paste(trait$Genus, trait$Species)
trait$FunGroup <- paste(trait$Native.Invasive, trait$Forb.Grass)

# 1995 is when all plots were done being added, but plots still changed sizes after this
# plot sizes were constant after 2001
sd <-  filter(sd, !species %in% remove, Year >= 2001, !seeds %in% c("-", "-99", "0"))
sd$seeds <- as.numeric(sd$seeds)

sd <- ddply(sd, .(Year, habitat, plot.habitat.replicate, species), summarize, PA = 1)
sd <- merge(sd, plot.size[,c(1,2,7)], by = c("Year", "plot.habitat.replicate"))

sd <- filter(sd, actual.plot.size == "0.1")

sd <- pivot_wider(sd[,1:5], names_from = "Year", values_from = "PA")
sd[is.na(sd)] <- 0
sd$'2006' <- 0 # 2006 was a complete failed germination year
sd$'2002' <- 0 # 2002 was a complete failed seed production year

length(unique(sd$plot.habitat.replicate)) #44

#### Create Species List ####
species.list <- list()

for(i in unique(sd$species)) {
  tmp <- filter(sd, species == i)
  tmp <- tmp[,-c(1:3)]
  tmp <- tmp[ , sort(colnames(tmp))]
  species.list[[i]] <- tmp
}

#### Run HMM ####

SD.df <- expand.grid(Species_Name = names(species.list), sim = 1:100, p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA, n.plots = NA, min.years = NA, max.years = NA, mean.years = NA) # empty df to fill with rates

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
    for (a in 1:5) print(i) 
    
    EMresult = EMestimation(X, r = 1) # update the model params, stop when new log likelihood is lower than the old log likelihood plus some precision 
   
    # fill in results with new params from the best log likelihood model
    p0Results[i] = EMresult$param$p0 # initial seed bank prob
    gResults[i] = EMresult$param$g # germ
    cResults[i] = EMresult$param$c #col
    sResults[i] = EMresult$param$s #surv
    rResults[i] = EMresult$param$r #prod
    llResults[i] = EMresult$ll # log-likelihood

  }
  
  SD.df[SD.df$Species_Name == j & SD.df$sim == k,]$p0 <- EMresult$param$p0
  SD.df[SD.df$Species_Name == j & SD.df$sim == k,]$g <- EMresult$param$g
  SD.df[SD.df$Species_Name == j & SD.df$sim == k,]$c <- EMresult$param$c
  SD.df[SD.df$Species_Name == j & SD.df$sim == k,]$s <- EMresult$param$s
  SD.df[SD.df$Species_Name == j & SD.df$sim == k,]$r <- EMresult$param$r
  SD.df[SD.df$Species_Name == j & SD.df$sim == k,]$iter <- which.max(EMresult$lllist)
  
  }

      cbind(j, 
            mean(SD.df[SD.df$Species_Name == j,]$p0), 
            mean(SD.df[SD.df$Species_Name == j,]$g), 
            mean(SD.df[SD.df$Species_Name == j,]$c), 
            mean(SD.df[SD.df$Species_Name == j,]$s),
            mean(SD.df[SD.df$Species_Name == j,]$r),
            mean(SD.df[SD.df$Species_Name == j,]$iter),
            sd(SD.df[SD.df$Species_Name == j,]$p0), 
            sd(SD.df[SD.df$Species_Name == j,]$g), 
            sd(SD.df[SD.df$Species_Name == j,]$c), 
            sd(SD.df[SD.df$Species_Name == j,]$s),
            sd(SD.df[SD.df$Species_Name == j,]$r),
            sd(SD.df[SD.df$Species_Name == j,]$iter,
            nrow(SD.df[which(SD.df$Species_Name == j & SD.df$iter < 150),])/100)
      )
}

#Stop cluster
stopCluster(cluster)

colnames(sim.df) <- c("Species_Name", "p0.boot", "g.boot", "c.boot", "s.boot", "r.boot", "iter.boot","p0.boot.sd", "g.boot.sd", "c.boot.sd", "s.boot.sd", "r.boot.sd", "iter.boot.sd", "perc.con")


saveRDS(sim.df, "SD-species-boot.RDS")
