# Correlation test for McLaughlin
rm(list=ls())

library(foreach)
library(doParallel)
library(dplyr)

source("~/Documents/00_HMM_Simulation.R")
source("~/Documents/00_HMM_Functions.R")

metadata <- read.csv("HMM-meta-mcl.csv")
mcl.df.q <- readRDS("McL-species-boot.RDS") # quadrat level HMM

mcl.df.q <- mcl.df.q %>% 
  dplyr::mutate(across(p0.boot:perc.con, ~as.numeric(.x)))

# update names for merger with trait data
mcl.df.q$Species_Name <- recode_factor(mcl.df.q$Species_Name, 'Lysimachia arvensis' = "Anagallis arvensis", 'Gastridium phleoides' = "Gastridium ventricosum")

mcl.df.q <- merge(mcl.df.q, metadata, by = "Species_Name")

mcl.df.q[!mcl.df.q$Species_Name %in% metadata$Species_Name,] # none

#### Filter data ####
mcl.df.q <- filter(mcl.df.q, n.plots >= 50) # 155 species to 62 species

mcl.df.q <- mcl.df.q[mcl.df.q$s.boot.sd < 0.1 & mcl.df.q$c.boot.sd < 0.1,] # 2 species lost

mcl.df.q <- mcl.df.q[mcl.df.q$perc.con > 0.5,] # zero species lost

mcl.df.q <- filter(mcl.df.q, FunGroup != "Native Grass") # one species lost

colnames(mcl.df.q)[2:6] <- c("p0", "g", "c", "s", "r")

#### simulation ###
param <- data.frame(matrix(ncol = 5, nrow = 5000))

colnames(param) <- c('p0', 'g', 'c', 's', 'r')

# independently sample from existing species parameters
param$p0 <- sample(mcl.df.q$p0, size = 5000, replace = T)
param$g <- sample(mcl.df.q$g, size = 5000, replace = T)
param$c <- sample(mcl.df.q$c, size = 5000, replace = T)
param$s <- sample(mcl.df.q$s, size = 5000, replace = T)
param$r <- 1

sim <- list()

for(i in 1:nrow(param)){
  sim[[i]] <- simulation(param[i,], N = 154, tf = 19)
}

#### Run HMM ####
trueParam = list() #object grouping together the parameters and other quantities which depend on them.
trueParam$c = 0.8      #colonization rate
trueParam$g = 0.5      #germination rate (sigma)
trueParam$r = 1        #reproductive rate (phi) keep this at 1 (if a seed germinates and survives to flower, we assume it produces a seed)
trueParam$s = 0.5      #seed survival
trueParam$p0 = 0.5     #initial state of the seed bank (probability that there were seeds in the soil the year before the first obs of existing flora)

trueParam = makeParametersCalculations(trueParam)

#Setup backend to use many processors
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-1) 
registerDoParallel(cluster)

sim.df <- foreach(j = 1:5000, .combine = rbind) %dopar% { # for each species
  X = as.matrix(sim[[j]]) # use their time series PA data
  n = 5 
  p0Results = rep(0,n)
  gResults = rep(0,n)
  cResults = rep(0,n)
  sResults = rep(0,n)
  rResults = rep(0,n)
  llResults = rep(0,n)
  lltrueParam = rep(0,n)
  
  for (i in 1:n) {
    EMresult = EMestimation(X, r = 1) # update the model params, stop when new log likelihood is lower than the old log likelihood plus some precision 
    
    # fill in results with new params from the best log likelihood model
    p0Results[i] = EMresult$param$p0 # initial seed bank prob
    gResults[i] = EMresult$param$g # germ
    cResults[i] = EMresult$param$c #col
    sResults[i] = EMresult$param$s #surv
    rResults[i] = EMresult$param$r #prod
    llResults[i] = EMresult$ll # log-likelihood
    lltrueParam[i] = logLikelihood(X,trueParam) # log likelihood of x given "true" parameters
  }
  
  cbind(j, 
        EMresult$param$p0, 
        EMresult$param$g, 
        EMresult$param$c, 
        EMresult$param$s,
        EMresult$param$r,
        which.max(EMresult$lllist))
}

#Stop cluster
stopCluster(cluster)

colnames(sim.df) <- c("Species_Name", "p0", "g", "c", "s", "r", "iter")

saveRDS(sim.df, "Intermediates/McL-sim-154-correlation.RDS")
