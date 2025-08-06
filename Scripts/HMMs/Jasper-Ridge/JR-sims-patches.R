# simulate Jasper Ridge
rm(list=ls())

library(foreach)
library(doParallel)
library(dplyr)

source("Scripts/HMMs/Validation/simulation.R")
source("Scripts/HMMs/HMM_Functions.R")

JR <- readRDS("Scripts/HMMs/Jasper-Ridge/JR-species-boot-avg.RDS")
JR <- as.data.frame(JR)

JR <- JR %>% 
  dplyr::mutate(across(p0.boot:r.boot, ~as.numeric(.x)))

# pre-filtered
# JR <- filter(JR, s.boot.sd < 0.1, c.boot.sd < 0.1, perc.con > 0.5)

rownames(JR) <- JR$Code

#### simulation ###
param <- JR[,1:6]
colnames(param) <- c("Species_Name", "p0", "g", "c", "s", "r")
rownames(param) <- param$Species_Name
param <- param[,-1]

#Setup backend to use many processors
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makeCluster(6) 
registerDoParallel(cluster)

tmp <- expand.grid(Species_Name = rownames(JR), sim = 1:100, p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA)

#### Run HMM ####
patch <- 15

sim.df <- foreach(j = rownames(JR), .combine = rbind, .errorhandling = "remove", .verbose = T) %dopar% { # for each species
  
  sim <- list()
  
  for(k in 1:5){
    sim[[k]] <- simulation(param[j,], N = patch, tf = 37)
    
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

saveRDS(sim.df, "SD-sim-species-validation-44.RDS")
