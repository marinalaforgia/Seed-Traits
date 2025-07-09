#### McLaughlin HMM ####
rm(list = ls())
source("~/Documents/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)


# Prep Data ####
trait <- read.csv("McL_80SitesSpeciesTraits_012615.csv")
mcl <- read.csv("Core_Community_Data2019.csv")
serp <- read.csv("McL_Abiotic-Data.csv")

# Confused and treated species (see data cleaning notes)
rm <- c("Logfia spp./Micropus californicus", "Aegilops triuncialis", "Trifolium bifidum/gracilentum", "Torilis arvensis/nodosa", "Hesperolinon spp.", "Cuscuta californica", "Silene gallica", "Erodium botrys", "Erodium brachycarpum", "Trifolium microdon", "Trifolium microcephalum") 

# rm commonly mis-ID'd species, species that hadnt always been censused (cuscuta) and goatgrass which is treated at the reserve and occasionally pulled from transects; 
# from susan email 8/9/22: "All Geranium should be dissectum"; "Silene gallica may be a mistake"
# from personal investigations: 
## EROBOT & EROBRA: Should ask susan but digging into the data there are plots where most are one, then they switch for one or two years, then go back (e.g., see site 36: EROBRA throughout but then EROBOT in 2010 and 2013). plus notes that say "should be botrys" but then aren't recorded as such
## Trifolium microdon & microcephalum: see sites 40, 59, 67 in 2013 for examples, notes saying ID uncertain/not flowering/could be microdon or microcephalum


mcl <- filter(mcl, Species_Name %in% trait[trait$Annual.Perennial == "Annual",]$Species_Name, !(Species_Name %in% rm)) 
mcl <- merge(serp[,1:2], mcl, by = "Site")

mcl[mcl$Species_Name == "Geranium molle",]$Species_Name <- "Geranium dissectum"
mcl <- ddply(mcl, .(Year, Site, Quadrat, Serpentine, Species_Name), summarize, Cover = sum(Cover))

##### Transect Level ####
# mcl <- ddply(mcl, .(Site, Year, Species_Name), summarize, PA = 1)
# mcl <- pivot_wider(mcl, names_from = "Year", values_from = "PA")
# mcl <- as.data.frame(mcl)
# mcl[is.na(mcl)] <- 0
# mcl <- merge(serp[,1:2], mcl, by = "Site")
# 
# mcl <- filter(mcl, Serpentine == "S")
# 
# # full <- data.frame()
# # for(s in names(sites)){
# #  mcl <- filter(mcl.tmp, Site %in% sites[[s]])
# 
# # looking at species averages
# mcl$avg.P.site <- apply(mcl[,3:22], 1, mean)
# mcl.sp.avg <- ddply(mcl, .(Species_Name), summarize, avg.P.sp = mean(avg.P.site))

##### Quadrat Level ####
#mcl[which(mcl$Cover == 0),] # checks out

mcl$PA <- 1
mcl$ID <- paste(mcl$Site, mcl$Quadrat, sep = ".")

mcl <- pivot_wider(mcl[,c(1,5,7,8)], names_from = "Year", values_from = "PA")
mcl[is.na(mcl)] <- 0

# looking at species averages
mcl$avg.P.site <- apply(mcl[,3:ncol(mcl)], 1, mean)
mcl.sp.avg <- ddply(mcl, .(Species_Name), summarize, avg.P.sp = mean(avg.P.site))

mcl <- mcl[,-23]

#### Create Species List ####
length(unique(mcl$ID))

#Full 400
# serp 190
# nonserp: 210
species.list <- list()

for(i in unique(mcl$Species_Name)) {
  tmp <- filter(mcl, Species_Name == i)
  tmp <- tmp[,-c(1:2)]
  tmp <- tmp[ , sort(colnames(tmp))]
  species.list[[i]] <- tmp 
}

#### Run HMM ####

mcl.df <- expand.grid(Species_Name = names(species.list), boot = 1:100, p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA, n.plots = NA, min.years = NA, max.years = NA, mean.years = NA) # empty df to fill with rates

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
  
  mcl.df[mcl.df$Species_Name == j & mcl.df$boot == k,]$p0 <- EMresult$param$p0
  mcl.df[mcl.df$Species_Name == j & mcl.df$boot == k,]$g <- EMresult$param$g
  mcl.df[mcl.df$Species_Name == j & mcl.df$boot == k,]$c <- EMresult$param$c
  mcl.df[mcl.df$Species_Name == j & mcl.df$boot == k,]$s <- EMresult$param$s
  mcl.df[mcl.df$Species_Name == j & mcl.df$boot == k,]$r <- EMresult$param$r
  mcl.df[mcl.df$Species_Name == j & mcl.df$boot == k,]$iter <- which.max(EMresult$lllist)
  
  }

  # here I assumed that the last estimated parameters, were what the model converged on but I'm not sure that's correct?
      cbind(j, 
            mean(mcl.df[mcl.df$Species_Name == j,]$p0), 
            mean(mcl.df[mcl.df$Species_Name == j,]$g), 
            mean(mcl.df[mcl.df$Species_Name == j,]$c), 
            mean(mcl.df[mcl.df$Species_Name == j,]$s),
            mean(mcl.df[mcl.df$Species_Name == j,]$r),
            mean(mcl.df[mcl.df$Species_Name == j,]$iter),
            sd(mcl.df[mcl.df$Species_Name == j,]$p0), 
            sd(mcl.df[mcl.df$Species_Name == j,]$g), 
            sd(mcl.df[mcl.df$Species_Name == j,]$c), 
            sd(mcl.df[mcl.df$Species_Name == j,]$s),
            sd(mcl.df[mcl.df$Species_Name == j,]$r),
            sd(mcl.df[mcl.df$Species_Name == j,]$iter),
            nrow(mcl.df[which(mcl.df$Species_Name == j & mcl.df$iter < 150),])/100
      )
}

#Stop cluster
stopCluster(cluster)

colnames(boot.df) <- c("Species_Name", "p0.boot", "g.boot", "c.boot", "s.boot", "r.boot", "iter.boot","p0.boot.sd", "g.boot.sd", "c.boot.sd", "s.boot.sd", "r.boot.sd", "iter.boot.sd", "perc.con")

saveRDS(boot.df, "McL-species-boot.RDS")
