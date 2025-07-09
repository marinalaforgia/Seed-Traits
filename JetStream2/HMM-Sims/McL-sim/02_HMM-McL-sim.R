# simulate mclaughlin to validate species estimates (Fig S2 prep)

rm(list=ls())

library(foreach)
library(doParallel)
library(dplyr)

source("Scripts/00_HMM_Simulation.R")
source("Scripts/00_HMM_Functions.R")

metadata <- read.csv("Data/HMM-meta-mcl.csv")
mcl <- readRDS("Data/McL-species-boot.RDS")
mcl <- as.data.frame(mcl)

mcl <- mcl %>% 
  dplyr::mutate(across(p0.boot:perc.con, ~as.numeric(.x)))

# update names for merger with trait data
mcl$Species_Name <- recode_factor(mcl$Species_Name, 'Lysimachia arvensis' = "Anagallis arvensis", 'Gastridium phleoides' = "Gastridium ventricosum")

mcl <- merge(mcl, metadata, by = "Species_Name")

mcl[!mcl$Species_Name %in% metadata$Species_Name,] # none

#### Filter data ####
mcl <- filter(mcl, n.plots >= 50) # 155 species to 62 species

mcl <- mcl[mcl$s.boot.sd < 0.1 & mcl$c.boot.sd < 0.1,] # 2 species lost

mcl <- mcl[mcl$perc.con > 0.5,] # zero species lost

mcl <- filter(mcl, FunGroup != "Native Grass") # one species lost

#### simulation ###
param <- mcl[,1:6]
colnames(param) <- c("Species_Name", "p0", "g", "c", "s", "r")
rownames(param) <- param$Species_Name
param <- param[,-1]

#Setup backend to use many processors
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-2) 
registerDoParallel(cluster)

tmp <- expand.grid(Species_Name = rownames(param), sim = 1:100, p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA)

#### Run HMM ####
patch <- 50

sim.df <- foreach(j = rownames(param), .combine = rbind, .errorhandling = "remove", .verbose = T) %dopar% { # for each species
  
  sim <- list()
  
  for(k in 1:100){
    sim[[k]] <- simulation(param[j,], N = patch, tf = 19)
    
    X = as.matrix(sim[[k]]) # use their time series PA data
    n = 5 
    p0Results = rep(0,n)
    gResults = rep(0,n)
    cResults = rep(0,n)
    sResults = rep(0,n)
    rResults = rep(0,n)
    llResults = rep(0,n)
    
    for (i in 1:n) {
      
      EMresult = EMestimation(X, r = 1) # update the model params, stop when new log likelihood is lower than the old log likelihood plus some precision 
      
      # fill in results with new params from the best log likelihood model
      p0Results[i] = EMresult$param$p0 # initial seed bank prob
      gResults[i] = EMresult$param$g # germ
      cResults[i] = EMresult$param$c #col
      sResults[i] = EMresult$param$s #surv
      rResults[i] = EMresult$param$r #prod
      llResults[i] = EMresult$ll # log-likelihood
      
    }
    
    tmp[tmp$Species_Name == j & tmp$sim == k,]$p0 <- EMresult$param$p0
    tmp[tmp$Species_Name == j & tmp$sim == k,]$g <- EMresult$param$g
    tmp[tmp$Species_Name == j & tmp$sim == k,]$c <- EMresult$param$c
    tmp[tmp$Species_Name == j & tmp$sim == k,]$s <- EMresult$param$s
    tmp[tmp$Species_Name == j & tmp$sim == k,]$r <- EMresult$param$r
    tmp[tmp$Species_Name == j & tmp$sim == k,]$iter <- which.max(EMresult$lllist)
    
  }
  
  cbind(patch, j, 
        tmp[tmp$Species_Name == j,]$sim,
        tmp[tmp$Species_Name == j,]$p0, 
        tmp[tmp$Species_Name == j,]$g, 
        tmp[tmp$Species_Name == j,]$c, 
        tmp[tmp$Species_Name == j,]$s, 
        tmp[tmp$Species_Name == j,]$r, 
        tmp[tmp$Species_Name == j,]$iter 
  )
}



#Stop cluster
stopCluster(cluster)

colnames(sim.df) <- c("patches", "Species_Name", "sim", "p0", "g", "c", "s", "r", "iter")

saveRDS(sim.df, "McL-sim-species-validation-50.RDS")
