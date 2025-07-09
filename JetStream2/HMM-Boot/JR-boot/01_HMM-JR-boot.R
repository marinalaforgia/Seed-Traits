#### Jasper Ridge HMM Boot ####
rm(list=ls())
source("~/Documents/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)

#### Prep Data ####
jr <- read.csv("JR_cover_updated.csv")[,-1] #weird first column

# normal filtering
unique(jr$species)
rm <- c("ROCK", "TRSP", "ARDO", "BARE", "POSC", "CHPO","STPU","MECA", "LOUT", "SIJU", "MISP", "THSP", "ESCA") # remove perennials and species that don't show up

jr[jr$species == "BRMO",]$species <- "BRHO"
jr[jr$species == "BR__",]$species <- "BRHO"
jr[jr$species == "BRTE",]$species <- "BRHO"
jr[jr$species == "BRTR",]$species <- "BRHO"
jr[jr$species == "EVSP",]$species <- "HESP"
jr[jr$species == "ORDE",]$species <- "CADE"
jr[jr$species == "EPPA",]$species <- "EPBR"
jr[jr$species == "HELU",]$species <- "HECO"
jr[jr$species == "LIAN",]$species <- "LIPA"
jr[jr$species == "LOSU",]$species <- "LOWR"
jr[jr$species == "TIER",]$species <- "CRCO"
jr[jr$species == "TRTR",]$species <- "TRWI" # doesnt work with 1m2

jr <-  filter(jr, !(species %in% rm))
jr <- jr %>%
  dplyr::group_by(quadID, species, year, treatment) %>%
  dplyr::summarize(cover = sum(cover)) # merge bromes in 0.5m plots

jr$PA <- ifelse(jr$cover > 0, 1, 0)

jr.sum <- pivot_wider(jr[,c(1:3,6)], names_from = "year", values_from = "PA")
jr.sum[is.na(jr.sum)] <- 0
jr.sum$sumPA <- rowSums(jr.sum[,3:39])
jr.sum <- filter(jr.sum, sumPA > 0) 

jr.sum <- filter(jr.sum, !species %in% rm)

control.plots <- c("c101", "c103", "c105", "c108", "c110", "c112", "c113", "c115", "c117", "c120", "c122", "c124", "c201", "c203", "c205", "c208", "c210", "c212", "c213", "c215", "c217", "c220", "c222", "c224", "c301", "c303", "c305", "c308", "c310", "c312", "c313", "c315", "c317", "c320", "c322", "c324")

#control.plots <- c("c102", "c104", "c106", "c107", "c109", "c111", "c114", "c116", "c118", "c119", "c121", "c123", "c202", "c204", "c206", "c207", "c209", "c211", "c214", "c216", "c218", "c219", "c221", "c223", "c302", "c304", "c306", "c307", "c309", "c311", "c314", "c316", "c318", "c319", "c321", "c323")

jr.sum <- filter(jr.sum, quadID %in% control.plots)
jr.sum <- jr.sum[,-40]

jr.sum <- filter(jr.sum, species != "PLER")

length(unique(jr.sum$quadID))

#### Create Species List ####
species.list <- list()

for(i in unique(jr.sum$species)) {
  tmp <- filter(jr.sum, species == i)
  tmp <- tmp[,-c(1:2)]
  tmp <- tmp[ , sort(colnames(tmp))]
  species.list[[i]] <- tmp 
}


#### Run HMM ####

JR.df <- expand.grid(Species_Name = names(species.list), boot = 1:100, p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA, n.plots = NA, min.years = NA, max.years = NA, mean.years = NA) # empty df to fill with rates

#Setup backend to use many processors
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-4) 
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
    
    JR.df[JR.df$Species_Name == j & JR.df$boot == k,]$p0 <- EMresult$param$p0
    JR.df[JR.df$Species_Name == j & JR.df$boot == k,]$g <- EMresult$param$g
    JR.df[JR.df$Species_Name == j & JR.df$boot == k,]$c <- EMresult$param$c
    JR.df[JR.df$Species_Name == j & JR.df$boot == k,]$s <- EMresult$param$s
    JR.df[JR.df$Species_Name == j & JR.df$boot == k,]$r <- EMresult$param$r
    JR.df[JR.df$Species_Name == j & JR.df$boot == k,]$iter <- which.max(EMresult$lllist)
    
  }
  
  cbind(j, 
        mean(JR.df[JR.df$Species_Name == j,]$p0), 
        mean(JR.df[JR.df$Species_Name == j,]$g), 
        mean(JR.df[JR.df$Species_Name == j,]$c), 
        mean(JR.df[JR.df$Species_Name == j,]$s),
        mean(JR.df[JR.df$Species_Name == j,]$r),
        mean(JR.df[JR.df$Species_Name == j,]$iter),
        sd(JR.df[JR.df$Species_Name == j,]$p0), 
        sd(JR.df[JR.df$Species_Name == j,]$g), 
        sd(JR.df[JR.df$Species_Name == j,]$c), 
        sd(JR.df[JR.df$Species_Name == j,]$s),
        sd(JR.df[JR.df$Species_Name == j,]$r),
        sd(JR.df[JR.df$Species_Name == j,]$iter),
        nrow(JR.df[which(JR.df$Species_Name == j & JR.df$iter < 150),])/100
  )
}

#Stop cluster
stopCluster(cluster)

colnames(boot.df) <- c("Species_Name", "p0.boot", "g.boot", "c.boot", "s.boot", "r.boot", "iter.boot","p0.boot.sd", "g.boot.sd", "c.boot.sd", "s.boot.sd", "r.boot.sd", "iter.boot.sd", "perc.con")

saveRDS(boot.df, "JR-species-boot-B.RDS")
