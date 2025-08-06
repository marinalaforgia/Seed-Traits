#### Portal HMM v2 ####
rm(list=ls())
source("Scripts/HMMs/HMM_Functions.R")
#use_default_data_path("Data/Long-Term-Datasets/Portal/")
#PORTALR_DATA_PATH="Data/Long-Term-Datasets/Portal/"

#### Load Libraries ####
library(ggplot2)
library(tidyverse)
library(portalr)

plant_data <- load_plant_data(
  path = "Data/Long-Term-Datasets/Portal/",
  download_if_missing = F,
  quiet = FALSE
)

plant.data <- clean_plant_data(
  plant_data,
  type = 'annual',
  unknowns = FALSE,
  correct_sp = TRUE
)


#plant.data <- filter(plant.data, season == "winter", community == "Winter Annual") 
#plant.data <- filter(plant.data, season == "summer", community == "Summer Annual")
plant.data <- filter(plant.data, community == "Summer and Winter Annual")

# only winter census and only winter SPECIES (no species that occur in both seasons as this will affect P/A but since we aren't looking at summer and they could be there in summer these species' estiamtes may be incorrect

# treatments change over time (See https://github.com/weecology/PortalData/blob/main/SiteandMethods/Portal_plot_treatments.csv) so although the project dates back to 1977, running models on controls must be split into different segments:

#controls <- c(4:7, 11, 13, 14, 17, 18, 24) #2015-2021 nvm apparently controls werent censused in 2020?
#controls <- c(1,2,9,11,14,22) #2005-2014 nvm unfortunately 2010 and 2013 not censused
controls <- c(2,11,14,22) #controls from 1988 - 2009; 2010 not censused
#no.kr <- c(1,5,6,7,9, 15, 16, 18, 21, 24)
#no.ant <- c(4, 8, 17, 12)
#no.kr.no.ant <- c(3, 10, 13, 19, 20, 23)

plant.data <- filter(plant.data, year >= 1988 & year < 2010, plot %in% controls)

#### Prep Data ####
portal <- plant.data

# overall mean
portal.RA <- portal %>%
  dplyr::mutate(sum.abundance = sum(abundance)) %>%
  dplyr::group_by(species, sum.abundance) %>%
  dplyr::summarize(sp.RA = sum(abundance)) %>%
  dplyr::group_by(species) %>%
  dplyr::summarize(RA = sp.RA/sum.abundance)

portal$PA <- 1

# for summer and winter species dont actually have two gens a year, they just flower in winter or later in summer

portal <- portal[,c(1,3:5,19)]
portal <- portal %>%
  dplyr::group_by(year, plot, quadrat, species) %>%
  dplyr::summarize(PA = mean(PA))

portal <- pivot_wider(portal, names_from = "year", values_from = "PA")

portal[is.na(portal)] <- 0

# according to the metadata 2009 were all done... https://github.com/weecology/PortalData/blob/main/Plants/Portal_plant_censuses.csv
portal$'2009' <- 0  # winter/summer

### for summer OR winter species
portal <- portal[,c(1,3:5,19)] # winter or summer only
portal <- portal %>%
  dplyr::group_by(year, plot, quadrat, species) %>%
  dplyr::summarize(PA = mean(PA))

portal <- pivot_wider(portal, names_from = "year", values_from = "PA")

portal[is.na(portal)] <- 0

# according to the metadata 1996, 2000, 2009 were all done... https://github.com/weecology/PortalData/blob/main/Plants/Portal_plant_censuses.csv
portal$'1996' <- 0   #winter ONLY
portal$'2000' <- 0  # winter ONLY
portal$'2009' <- 0  # winter/summer
portal$'2003' <- 0  # SUMMER ONLY

####


portal <- portal[ , order(names(portal))]
portal <- portal[,c(23:25, 1:22)] 
colnames(portal)[3] <- "SpeciesCode"

#### Metadata ####
portal$sumPA <- rowSums(portal[,4:ncol(portal)])

meta.port <- portal %>%
  dplyr::group_by(SpeciesCode) %>%
  dplyr::summarize(n.plots = length(plot),
            mean.years = mean(sumPA),
            min.years = min(sumPA),
            max.years = max(sumPA))

portal <- portal[,-ncol(portal)]

#### Create Species List ####
length(unique(paste(portal$plot, portal$quadrat))) # 64 control plots, 160 no kangaroo rat plots, 48 no ant, 96 no ant no rod; 384 combined

species.list <- list()

for(i in unique(portal$SpeciesCode)) {
  tmp <- filter(portal, SpeciesCode == i)
  tmp <- tmp[,-c(1:3)]
  tmp <- tmp[ , sort(colnames(tmp))]
  species.list[[i]] <- tmp 
}


#### Run HMM ####

port.df <- expand.grid(Species_Name = names(species.list), p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA, n.plots = NA, min.years = NA, max.years = NA, mean.years = NA) # empty df to fill with rates

for(j in names(species.list)){ # for each species
  X = as.matrix(species.list[[j]]) # use their time series PA data
  n = 5 # why n = 5?
  p0Results = rep(0,n)
  gResults = rep(0,n)
  cResults = rep(0,n)
  sResults = rep(0,n)
  rResults = rep(0,n)
  llResults = rep(0,n)
  
  for (i in 1:n) {
    for (k in 1:5) print(i) 
    EMresult = EMestimation(X, r = 1) # update the model params, stop when new log likelihood is lower than the old log likelihood plus some precision 
  
    # fill in results with new params from the best log likelihood model
    p0Results[i] = EMresult$param$p0 # initial seed bank prob
    gResults[i] = EMresult$param$g # germ
    cResults[i] = EMresult$param$c #col
    sResults[i] = EMresult$param$s #surv
    rResults[i] = EMresult$param$r #prod
    llResults[i] = EMresult$ll # log-likelihood

  }

  # last estimated parameters ie what the model converged on
  port.df[port.df$Species_Name == j,]$p0 <- EMresult$param$p0
  port.df[port.df$Species_Name == j,]$g <- EMresult$param$g
  port.df[port.df$Species_Name == j,]$c <- EMresult$param$c
  port.df[port.df$Species_Name == j,]$s <- EMresult$param$s
  port.df[port.df$Species_Name == j,]$r <- EMresult$param$r
  port.df[port.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
  port.df[port.df$Species_Name == j,]$n.plots <- nrow(species.list[[j]])
  port.df[port.df$Species_Name == j,]$mean.years <- mean(rowSums(species.list[[j]]))
  port.df[port.df$Species_Name == j,]$min.years <- min(rowSums(species.list[[j]]))
    port.df[port.df$Species_Name == j,]$max.years <- max(rowSums(species.list[[j]]))
  
}

#### Explore Output ####
trait.repo <- read.csv("Data/Long-Term-Datasets/Portal/repo-traits.csv")
trait.bio <- read.csv("Data/Long-Term-Datasets/Portal/biotime-portal-species.csv")

traits <- merge(trait.repo, trait.bio, by = "Code", all = T)
traits$Species_Name <- paste(traits$genus, traits$sp, sep = " ")

colnames(port.df)[1] <- "species"
port.df <- merge(port.df, traits[, c(1,2,13,17,18)], by = "species", all.x = T, all.y = F)
port.df <- filter(port.df, Code != "ASTALL", Code != "BAIMUL") # remove  perennials ASTALL and BAIMUL


#port.df[port.df$Species_Name == "Dimorphocarpa wislizeni",]$FunGroup <- "Native Forb"

port.df[port.df$Code == "CHEFRE",]$Code <- "CHEINC"

port.df$FunGroup <- paste(port.df$Native.Invasive, port.df$Grass.Forb, sep = " ")

port.df <- port.df[,c(13,12,16,2:11)]

ggplot(port.df[port.df$iter<150,], aes(x = s, y = c)) +
  #geom_point() +
  geom_text(aes(label = Code)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_classic()

#saveRDS(port.df, "Scripts/HMMs/Portal/Port_HMM_control_winter_1988-2009.RDS")
#saveRDS(port.df, "Scripts/HMMs/Portal/Port_HMM_control_summer_1988-2009.RDS")
saveRDS(port.df, "Scripts/HMMs/Portal/Port_HMM_control_summer-winter_1988-2009.RDS")

#### change names in RA data ####
portal.RA <- merge(portal.RA, traits[, c(2,1)], by = "species", all.x = T, all.y = F)
portal.RA <- portal.RA[,c(3,2)]
saveRDS(portal.RA, "Data/Long-Term-Datasets/Portal/Portal_RA.RDS")

#### Metadata ####
meta.port <- merge(meta.port, traits[, c(1,2,13,17,18)], by.x = "SpeciesCode", by.y = "species", all.x = T, all.y = F)
meta.port$FunGroup <- paste(meta.port$Native.Invasive, meta.port$Grass.Forb, sep = " ")

meta.port <- filter(meta.port, Code != "ASTALL", Code != "BAIMUL")
meta.port <- meta.port[,c(1,7,6,10,2,4,5,3)]
write.csv(meta.port, "Data/Long-Term-Datasets/Portal/HMM-meta-port-summer-winter.csv", row.names = F)

port.df.s <- readRDS("Scripts/HMMs/Portal/Archive/Port_HMM_control_summer_1988-2009.RDS")
port.df.w <- readRDS("Scripts/HMMs/Portal/Archive/Port_HMM_control_winter_1988-2009.RDS")
port.df.sw <- readRDS("Scripts/HMMs/Portal/Archive/Port_HMM_control_summer-winter_1988-2009.RDS")

port.df <- rbind(port.df.s, port.df.w, port.df.sw)
port.df[port.df$Code == "ERAPEC",]$FunGroup <- "Native Grass"
port.df[port.df$Code == "CHANIC",]$FunGroup <- "Native Forb"

ggplot(port.df[port.df$iter < 150 & port.df$n.plots > 10,], aes(x = s, y = c)) +
  #geom_point() +
  geom_text(aes(label = Code)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_classic()
