#### McLaughlin HMM ####
rm(list = ls())
source("Scripts/HMMs/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(plyr)
library(tidyverse)

# Prep Data ####
skin <- read.csv("Data/Long-Term-Datasets/Skinner/L Skinner unburned 1994-2005.csv")
sp <- read.csv("Data/Long-Term-Datasets/Skinner/Species-List.csv")
rm <- c("Bareground", "LITTER","Erodium spp","Galium spp.","Artemisia californica sdlg","Encelia farinosa seedling","Eriodyction crassifolium","Eriogonum fasciculatum sdlg","Salvia mellifera sdlg.","Sums","Sum EF","Sum EG","Sum NF","Sum NG","Sum NS","Calochortus venustus?","Camissonia sp.","Hemizonia sp.","Navarretia sp", "Dichelostemma capitata", "Stipa pulchra", "Acmispon glaber", "Brassica geniculata") # rm unknowns, perennials

skin <- pivot_longer(skin, cols = Brassica.geniculata:Sum.NS, names_to = "Species.Name", values_to = "Cover")

skin <- merge(skin, sp[,c(1,3,4)], by = "Species.Name")

skin <- filter(skin, Cover > 0, !Species_Name %in% rm)

#skin <- filter(skin, F.NF == "NF")
skin <- skin[,c(2,7,12,10)]
colnames(skin)[c(1,2)] <- c("Year", "Quadrat")
skin <- ddply(skin, .(Year, Quadrat, Code), summarize, Cover = sum(Cover))
skin$PA <- 1

skin.RA <- skin %>%
  dplyr::mutate(sum.abundance = sum(Cover)) %>%
  dplyr::group_by(Code, sum.abundance) %>%
  dplyr::summarize(sp.RA = sum(Cover)) %>%
  dplyr::group_by(Code) %>%
  dplyr::summarize(RA = sp.RA/sum.abundance)

saveRDS(skin.RA, "Scripts/HMMs/Skinner/skin-ra.RDS")

skin <- pivot_wider(skin[,c(1:3,5)], names_from = "Year", values_from = "PA")
skin[is.na(skin)] <- 0
skin$'2002' <- 0

#### Create Species List ####
length(unique(skin$Quadrat)) #16

species.list <- list()
# annoyingly sensitive to the first year of data - cant be all 1s for a species or the model breaks, if it's all zeroes, s is estimated at 0
skin <- filter(skin, Code != "BRORUB")
               
for(i in unique(skin$Code)) {
  tmp <- filter(skin, Code == i)
  tmp <- tmp[,-c(1:2)]
  tmp <- tmp[ , sort(colnames(tmp))]
  if(nrow(tmp) > 10) species.list[[i]] <- tmp  # for now I got rid of species that were present in less than 10 sites throughout the entire time series
}


df <- data.frame(names = NA, obs = NA)

for(n in names(species.list)) {
  df <- data.frame(rbind(df, c(names = n, obs = nrow(species.list[[n]]))))
}


df$obs <- as.numeric(as.character(df$obs))


#### Run HMM ####
# starting values. might not necessarily help much to change these, but we can try
trueParam = list() #object grouping together the parameters and other quantities which depend on them.
trueParam$c = 0.8      #colonization rate
trueParam$g = 0.5      #germination rate (sigma)
trueParam$r = 1        #reproductive rate (phi) keep this at 1 (if a seed germinates and survives to flower, we assume it produces a seed)
trueParam$s = 0.5      #seed survival
trueParam$p0 = 0.5     #initial state of the seed bank (probability that there were seeds in the soil the year before the first obs of existing flora)

trueParam = makeParametersCalculations(trueParam)


skin.df <- expand.grid(Species_Name = names(species.list), p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA, n.plots = NA, min.years = NA, max.years = NA, mean.years = NA) # empty df to fill with rates

# Took about 20 minutes or so
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
 
  print(j)  
  
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
  skin.df[skin.df$Species_Name == j,]$p0 <- EMresult$param$p0
  skin.df[skin.df$Species_Name == j,]$g <- EMresult$param$g
  skin.df[skin.df$Species_Name == j,]$c <- EMresult$param$c
  skin.df[skin.df$Species_Name == j,]$s <- EMresult$param$s
  skin.df[skin.df$Species_Name == j,]$r <- EMresult$param$r
  skin.df[skin.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
  skin.df[skin.df$Species_Name == j,]$n.plots <- nrow(species.list[[j]])
  skin.df[skin.df$Species_Name == j,]$mean.years <- mean(rowSums(species.list[[j]]))
  skin.df[skin.df$Species_Name == j,]$min.years <- min(rowSums(species.list[[j]]))
  skin.df[skin.df$Species_Name == j,]$max.years <- max(rowSums(species.list[[j]]))
}

traits <- read.csv("Data/20211001_Full-Species-List.csv")
traits$FunGroup <- paste(traits$nat.inv, traits$group, sep = " ")

colnames(skin.df)[1] <- "Code"
skin.df <- merge(skin.df, sp[,c(3,4)], by = "Code", all.x = T, all.y = F)

skin.df <- merge(skin.df, traits[,c(7,22)], by.x = "Code", by.y = "new.code", all.x = T, all.y = F)

skin.df[skin.df$Code == "ACMMIC",]$FunGroup <- "native forb"
skin.df[skin.df$Code == "CAMBIS",]$FunGroup <- "native forb"
skin.df[skin.df$Code == "EUPALB",]$FunGroup <- "native forb"
skin.df[skin.df$Code == "MIMBRE",]$FunGroup <- "native forb"
skin.df[skin.df$Code == "PECLIN",]$FunGroup <- "native forb"
skin.df[skin.df$Code == "STYGNA",]$FunGroup <- "native forb"
skin.df[skin.df$Code == "CRACON",]$FunGroup <- "native forb"

ggplot(skin.df[skin.df$iter < 150,], aes(x = s, y = c)) +
  geom_point() +
  geom_text(aes(label = Code)) + 
  geom_smooth(method = "lm", formula = y ~ x, se = F) +
  theme_bw()


skin.df$FunGroup <- ifelse(skin.df$FunGroup == "native forb", "Native Forb", ifelse(skin.df$FunGroup == "native grass", "Native Grass", skin.df$FunGroup))

skin.df$FunGroup <- ifelse(skin.df$FunGroup == "invasive forb", "Exotic Forb", ifelse(skin.df$FunGroup == "invasive grass", "Exotic Grass", skin.df$FunGroup))

skin.df <- skin.df[,c(12,1,13,2:11)]

saveRDS(skin.df, "Scripts/HMMs/Skinner/skinner-hmm.RDS")



