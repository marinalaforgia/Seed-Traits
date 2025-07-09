# simulate mclaughin
rm(list=ls())

library(foreach)
library(doParallel)
library(plyr)
library(dplyr)

source("simulation.R")
source("HMM_Functions.R")

mcl.df.q <- readRDS("20221212_mcl_HMM_quadrat.RDS") # quadrat level HMM
mcl.df.q <- filter(mcl.df.q, iter < 150, n.plots >= 50, FunGroup != "Native Grass")

#### simulation ###
param <- data.frame(matrix(ncol = 5, nrow = 5000))

colnames(param) <- c('p0', 'g', 'c', 's', 'r')

param$p0 <- sample(mcl.df.q$p0, size = 5000, replace = T)
param$g <- sample(mcl.df.q$g, size = 5000, replace = T)
param$c <- sample(mcl.df.q$c, size = 5000, replace = T)
param$s <- sample(mcl.df.q$s, size = 5000, replace = T)
param$r <- 1

sim <- list()

for(i in 1:nrow(param)){
  sim[[i]] <- simulation(param[i,])
}

#### Run HMM ####


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

saveRDS(sim.df, "McL-sim-154.RDS")
