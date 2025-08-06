#### Portal HMM ####
rm(list=ls())

source("HMMs/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(tidyverse)

#### Winter ####
#### .___Prep Data ####
portal <- read.csv("Data/Long-Term-Datasets/Portal/Portal_plant_winter_annual_19892002.csv")
trt <- read.csv("Data/Long-Term-Datasets/Portal/Treatment.csv")
rm <- c("D_sp2", "Cryp_sp3", "Cryp_sp4", "Delp_sp1", "D_sp2", "Lupi_sp", "Lupi_sp2", "Oeno_sp", "Pect_sp", "unk1", "unk2", "unk3", "unk4", "unk5", "Desc_sp2", "ukn_grass","Astr_spp") # remove unknowns

portal <- portal[,!(colnames(portal) %in% rm)]
#portal <- filter(portal, Year %in% 1989:1997) # date range used in Chen & Valone 2017 as after this Erodium takes over... the models did not seem to like this however...
portal <- portal %>% mutate(across(6:49, ~1 * (. > 0)))

portal <- merge(portal, trt[c(1,5)], by = "Plot")

portal <- filter(portal, Trt == "C")

portal$plot.quad <- paste(portal$Plot, portal$Quadrat, sep = ".")
portal <- portal[,c(51,4,6:49)] #with treat

portal <- portal %>% 
  mutate(across(Ambr_arte:Vulp_octo, as.numeric)) %>% 
  pivot_longer(names_to = "SpeciesCode", cols = Ambr_arte:Vulp_octo) %>% 
  pivot_wider(id_cols = c(plot.quad, SpeciesCode), values_from = value,  names_from = Year)

portal$sum <- rowSums(portal[,3:11])
portal <- filter(portal, sum > 0)
portal <- portal[,-12]

portal <- portal[ , order(names(portal))]
portal <- portal[,c(10,11,1:9)]

#### .___Create Species List ####

species.list <- list()

for(i in unique(portal$SpeciesCode)) {
  tmp <- filter(portal, SpeciesCode == i)
  tmp <- tmp[,-c(1:2)]
  tmp <- tmp[ , sort(colnames(tmp))]
  if(nrow(tmp) > 8) species.list[[i]] <- tmp 
}

df <- data.frame()

for(i in names(species.list)) {
  tmp <- data.frame(names = i, obs = nrow(species.list[[i]]))
  df <- rbind(df, tmp)
}

df$obs <- as.numeric(as.character(df$obs))


#### .___Run HMM ####
# starting values. might not necessarily help much to change these, but we can try
trueParam = list() #object grouping together the parameters and other quantities which depend on them.
trueParam$c = 0.2      #colonization rate
trueParam$g = 0.5      #germination rate (sigma)
trueParam$r = 1        #reproductive rate (phi) keep this at 1 (if a seed germinates and survives to flower, we assume it produces a seed)
trueParam$s = 0.5      #seed survival
trueParam$p0 = 0.5     #initial state of the seed bank (probability that there were seeds in the soil the year before the first obs of existing flora)

trueParam = makeParametersCalculations(trueParam)


port.df <- expand.grid(Species_Name = names(species.list), p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA) # empty df to fill with rates

for(j in names(species.list)){ # for each species
  X = as.matrix(species.list[[j]]) # use their time series PA data
  n = 5 # why n = 5?
  p0Results = rep(0,n)
  gResults = rep(0,n)
  cResults = rep(0,n)
  sResults = rep(0,n)
  rResults = rep(0,n)
  llResults = rep(0,n)
  lltrueParam = rep(0,n)
  
  for (i in 1:n) {
    for (k in 1:5) print(i) # I don't understand what's going on here, why print it 5 times?
    print(logLikelihood(X, trueParam)) # print the log-likelihood given the "true params"
    EMresult = EMestimation(X, r = 1) # update the model params, stop when new log likelihood is lower than the old log likelihood plus some precision 
    print(logLikelihood(X, trueParam)) # print the log-likelihood given the "true params"; why is this on here twice?
    # fill in results with new params from the best log likelihood model
    p0Results[i] = EMresult$param$p0 # initial seed bank prob
    gResults[i] = EMresult$param$g # germ
    cResults[i] = EMresult$param$c #col
    sResults[i] = EMresult$param$s #surv
    rResults[i] = EMresult$param$r #prod
    llResults[i] = EMresult$ll # log-likelihood
    lltrueParam[i] = logLikelihood(X,trueParam) # log likelihood of x given "true" parameters
  }

  # here I assumed that the last estimated parameters, were what the model converged on but I'm not sure that's correct?
  port.df[port.df$Species_Name == j,]$p0 <- EMresult$param$p0
  port.df[port.df$Species_Name == j,]$g <- EMresult$param$g
  port.df[port.df$Species_Name == j,]$c <- EMresult$param$c
  port.df[port.df$Species_Name == j,]$s <- EMresult$param$s
  port.df[port.df$Species_Name == j,]$r <- EMresult$param$r
  port.df[port.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
  
}

# remove species that didn't converge in 150 steps
#port.df <- filter(port.df, iter < 100)

#### .___Explore Output ####
trait <- read.csv("Data/Long-Term-Datasets/Portal/biotime-portal-species.csv")
port.df <- merge(port.df, trait[, c(2,3,5,6,7)], by.x = "Species_Name", by.y = "SpeciesCode", all.x = T, all.y = F)
port.df <- filter(port.df, include == "y")
port.df$FunGroup <- paste(port.df$Native.Invasive, port.df$Grass.Forb, sep = " ")

colnames(port.df)[1] <- "Code"
colnames(port.df)[8] <- "Species_Name"
# #port.df <- port.df[,c(9,1,13,2:7,8)]
port.df <- port.df[,c(8,1,12,2:7)]

ggplot(port.df[port.df$iter<100,], aes(x = s, y = c)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_classic()

saveRDS(port.df, "HMMs/Portal/Port_winter_all-controls_HMM.RDS")

#### Summer  ####
#### .___Prep Data ####
portal <- read.csv("Data/Long-Term-Datasets/Portal/Portal_plant_summer_annual_19892002.csv")
trt <- read.csv("Data/Long-Term-Datasets/Portal/Treatment.csv")
rm <- c("Bout_sp", "Boer_sp", "Cusc_sp", "Pani_sp", "grass_unk", "Port_sp", "Sida_sp", "unkn_sp", "Moll_sp") # remove unknowns

portal <- portal[,!(colnames(portal) %in% rm)]

portal <- portal %>% mutate(across(6:54, ~1 * (. > 0)))

portal <- merge(portal, trt[c(1,5)], by = "Plot")

portal <- filter(portal, Trt == "C")

portal <- filter(portal, Quadrat %in% c(11,33,55,77,15,37,51,73))

portal$plot.quad <- paste(portal$Plot, portal$Quadrat, sep = ".")
portal <- portal[,c(56,4,6:54)] #with treat

portal <- portal %>% 
  mutate(across(Amar_palm:Verb_ence, as.numeric)) %>% 
  pivot_longer(names_to = "SpeciesCode", cols = Amar_palm:Verb_ence) %>% 
  pivot_wider(id_cols = c(plot.quad, SpeciesCode), values_from = value,  names_from = Year)


portal$sum <- rowSums(portal[,3:16])
portal <- filter(portal, sum > 0)
portal <- portal[,-17]

portal <- portal[ , order(names(portal))]
portal <- portal[,c(15,16,1:14)]

#### .___Create Species List ####

species.list <- list()

for(i in unique(portal$SpeciesCode)) {
  tmp <- filter(portal, SpeciesCode == i)
  tmp <- tmp[,-c(1:2)]
  tmp <- tmp[ , sort(colnames(tmp))]
  if(nrow(tmp) > 8) species.list[[i]] <- tmp 
}

df <- data.frame()

for(i in names(species.list)) {
  tmp <- data.frame(names = i, obs = nrow(species.list[[i]]))
  df <- rbind(df, tmp)
}

df$obs <- as.numeric(as.character(df$obs))


#### .___Run HMM ####
# starting values. might not necessarily help much to change these, but we can try
trueParam = numeric(0) #object grouping together the parameters and other quantities which depend on them.
trueParam$c = 0.2      #colonization rate
trueParam$g = 0.5      #germination rate (sigma)
trueParam$r = 1        #reproductive rate (phi) keep this at 1 (if a seed germinates and survives to flower, we assume it produces a seed)
trueParam$s = 0.5      #seed survival
trueParam$p0 = 0.5     #initial state of the seed bank (probability that there were seeds in the soil the year before the first obs of existing flora)

trueParam = makeParametersCalculations(trueParam)


port.df <- expand.grid(Species_Name = names(species.list), p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA) # empty df to fill with rates

for(j in names(species.list)){ # for each species
  X = as.matrix(species.list[[j]]) # use their time series PA data
  n = 5 # why n = 5?
  p0Results = rep(0,n)
  gResults = rep(0,n)
  cResults = rep(0,n)
  sResults = rep(0,n)
  rResults = rep(0,n)
  llResults = rep(0,n)
  lltrueParam = rep(0,n)
  
  for (i in 1:n) {
    for (k in 1:5) print(i) # I don't understand what's going on here, why print it 5 times?
    print(logLikelihood(X, trueParam)) # print the log-likelihood given the "true params"
    EMresult = EMestimation(X, r = 1) # update the model params, stop when new log likelihood is lower than the old log likelihood plus some precision 
    print(logLikelihood(X, trueParam)) # print the log-likelihood given the "true params"; why is this on here twice?
    # fill in results with new params from the best log likelihood model
    p0Results[i] = EMresult$param$p0 # initial seed bank prob
    gResults[i] = EMresult$param$g # germ
    cResults[i] = EMresult$param$c #col
    sResults[i] = EMresult$param$s #surv
    rResults[i] = EMresult$param$r #prod
    llResults[i] = EMresult$ll # log-likelihood
    lltrueParam[i] = logLikelihood(X,trueParam) # log likelihood of x given "true" parameters
  }

  # here I assumed that the last estimated parameters, were what the model converged on but I'm not sure that's correct?
  port.df[port.df$Species_Name == j,]$p0 <- EMresult$param$p0
  port.df[port.df$Species_Name == j,]$g <- EMresult$param$g
  port.df[port.df$Species_Name == j,]$c <- EMresult$param$c
  port.df[port.df$Species_Name == j,]$s <- EMresult$param$s
  port.df[port.df$Species_Name == j,]$r <- EMresult$param$r
  port.df[port.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
  
}

# remove species that didn't converge in 150 steps
#port.df <- filter(port.df, iter < 100)

#### .___Explore Output ####
trait <- read.csv("Data/Long-Term-Datasets/Portal/biotime-portal-species.csv")
port.df <- merge(port.df, trait[, c(2,3,5,6,7)], by.x = "Species_Name", by.y = "SpeciesCode", all.x = T, all.y = F)
port.df <- filter(port.df, include == "y")
port.df$FunGroup <- paste(port.df$Native.Invasive, port.df$Grass.Forb, sep = " ")

colnames(port.df)[1] <- "Code"
colnames(port.df)[8] <- "Species_Name"
#port.df <- port.df[,c(9,1,13,2:7,8)]
port.df <- port.df[,c(8,1,12,2:7)]

ggplot(port.df[port.df$iter<100,], aes(x = s, y = c)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_classic()

saveRDS(port.df, "HMMs/Portal/Port_summer_sub-controls_HMM.RDS")
