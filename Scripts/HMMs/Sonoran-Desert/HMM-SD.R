#### Sonoran Desert HMM ####
rm(list=ls())
source("Scripts/HMMs/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(plyr)
library(tidyverse)

#### Prep Data ####
sd <- read.csv("Data/Long-Term-Datasets/Sonoran-Desert/CensusData.csv")
#m2 <- read.csv("Data/Long-Term-Datasets/Sonoran-Desert/LxBxperPlot_12sppREduced.csv")
#m2 <- unique(m2[,c(1,3,10)])

# Calculate plot sizes
# sd$PA <- 1
# sd.plot <- ddply(sd, .(Year, habitat, plot.habitat.replicate), summarize, seedlings = sum(PA))
# 
# sd.plot <- merge(sd.plot, m2, by.x = c("Year" ,"plot.habitat.replicate"), by.y = c("Year", "PlotID"), all.x = T)
# 
# sd.plot$plot.size <- round(sd.plot$seedlings/sd.plot$AllSppSPM2, 3)

# Larry does not have this information anywhere, but looking at counts over time these look to me like the 11 extra sites that weren't censused over time (this leaves the full 72 sites that were); Jenny confirmed this was the case
iffy.sites <- c("3-open-west", "3-open-east", "0/3-shrub-north", "30-open-north", "30-open-south", "30-shrub-north", "30-shrub-south", "7-open-east", "7-open-west", "8-open-east", "8-open-west")

#remove weird sites
sd <- filter(sd, !(plot.habitat.replicate %in% iffy.sites))
length(unique(sd$plot.habitat.replicate))
# write.csv(sd.plot, "Data/Long-Term-Datasets/Sonoran-Desert/plot-sizes.csv", row.names = F)
plot.size <- read.csv("Data/Long-Term-Datasets/Sonoran-Desert/plot-sizes.csv")


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



trait <- read.csv("Data/Long-Term-Datasets/Sonoran-Desert/Species_list_2016_2.csv")
trait$Species_Name <- paste(trait$Genus, trait$Species)
trait$FunGroup <- paste(trait$Native.Invasive, trait$Forb.Grass)


# 1995 is when all plots were done being added, but plots still changed sizes after this
# sd <-  filter(sd, !species %in% remove, Year >= 1995, !seeds %in% c("-", "-99", "0"))
# sd$seeds <- as.numeric(sd$seeds)

#plot sizes were constant after 2001
sd <-  filter(sd, !species %in% remove, Year >= 2001, !seeds %in% c("-", "-99", "0"))
sd$seeds <- as.numeric(sd$seeds)

# currently no way just looking at presence absence to adjust for plot size over time!  so we either have to assume if it wasn't present in 0.05 that it wouldnt be present in 0.1, or vice versa
sd <- ddply(sd, .(Year, habitat, plot.habitat.replicate, species), summarize, PA = 1)
sd <- merge(sd, plot.size[,c(1,2,7)], by = c("Year", "plot.habitat.replicate"))

sd <- filter(sd, actual.plot.size == "0.1")

length(unique(sd$plot.habitat.replicate))

sd.RA <- sd %>%
  dplyr::mutate(sum.abundance = sum(PA)) %>%
  dplyr::group_by(species, sum.abundance) %>%
  dplyr::summarize(sp.RA = sum(PA)) %>%
  dplyr::group_by(species) %>%
  dplyr::summarize(RA = sp.RA/sum.abundance)

sd.RA <- merge(sd.RA, trait[,3:4], by.x = "species", by.y = "Species.code", all.x = T, all.y = F)
sd.RA <- sd.RA[,c(3,2)]
saveRDS(sd.RA, "Data/Long-Term-Datasets/Sonoran-Desert/SD_RA_size01.RDS")

sd <- pivot_wider(sd[,1:5], names_from = "Year", values_from = "PA")
sd[is.na(sd)] <- 0
sd$'2006' <- 0 # 2006 was a complete failed germination year
sd$'2002' <- 0 # 2002 was a complete failed seed production year

length(unique(sd$plot.habitat.replicate)) #44

#### Metadata ####
sd$sumPA <- rowSums(sd[,4:19])

meta.sd <- sd %>%
  dplyr::group_by(species) %>%
  dplyr::summarize(n.plots = length(plot.habitat.replicate),
            mean.years = mean(sumPA),
            min.years = min(sumPA),
            max.years = max(sumPA))

sd <- sd[,-20]

meta.sd <- merge(meta.sd, trait[,c(3,4,13,14)], by.x = "species", by.y = "Species.code")

meta.sd <- meta.sd[,c(1,7,6,8,2,4,5,3)]

write.csv(meta.sd, "Data/Long-Term-Datasets/Sonoran-Desert/HMM-meta-SD.csv", row.names = F)

#### Create Species List ####
species.list <- list()

for(i in unique(sd$species)) {
  tmp <- filter(sd, species == i)
  tmp <- tmp[,-c(1:3)]
  tmp <- tmp[ , sort(colnames(tmp))]
  species.list[[i]] <- tmp
}

#### Run HMM ####

sd.df <- expand.grid(Species_Name = names(species.list), p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA, n.plots = NA, min.years = NA, max.years = NA, mean.years = NA) # empty df to fill with rates (iter = how many iterations it took the model to converge)

for(j in names(species.list)){ # for each species
  X = as.matrix(species.list[[j]]) # use their time series PA data
  n = 5 # why n = 5?
  p0Results = rep(0,n)
  gResults = rep(0,n)
  cResults = rep(0,n)
  sResults = rep(0,n)
  rResults = rep(0,n)
  llResults = rep(0,n)

  print(j)  
 
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

  sd.df[sd.df$Species_Name == j,]$p0 <- EMresult$param$p0
  sd.df[sd.df$Species_Name == j,]$g <- EMresult$param$g
  sd.df[sd.df$Species_Name == j,]$c <- EMresult$param$c
  sd.df[sd.df$Species_Name == j,]$s <- EMresult$param$s
  sd.df[sd.df$Species_Name == j,]$r <- EMresult$param$r
  sd.df[sd.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
  sd.df[sd.df$Species_Name == j,]$n.plots <- nrow(species.list[[j]])
  sd.df[sd.df$Species_Name == j,]$mean.years <- mean(rowSums(species.list[[j]]))
  sd.df[sd.df$Species_Name == j,]$min.years <- min(rowSums(species.list[[j]]))
  sd.df[sd.df$Species_Name == j,]$max.years <- max(rowSums(species.list[[j]]))
}

#### Explore Output ####
colnames(sd.df)[1] <- "code.old"
sd.df <- merge(sd.df, trait[,c(3,4,13,14)], by.x = "code.old", by.y = "Species.code")

sd.df <- sd.df[,c(13,12,14,2:11)]

saveRDS(sd.df, "Scripts/HMMs/Sonoran-Desert/SD_HMM_2001_01.RDS")

ggplot(sd.df[sd.df$iter < 150,], aes(x = s, y = c)) +
  #geom_point() +
  geom_text(aes(label = Code)) +
  theme_classic() +
  geom_smooth(method = "lm", formula = y ~ x)

