#### McLaughlin HMM ####
rm(list = ls())
source("Scripts/HMMs/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(plyr)
library(tidyverse)

# Prep Data ####
trait <- read.csv("/Users/Marina/Google Drive/02_McLaughlin_80Sites_Organized/Modified_CombinedFiles/McL_80SitesSpeciesTraits_012615.csv")
mcl <- read.csv("/Users/Marina/Google Drive/02_McLaughlin_80Sites_Organized/Data-Storage-Files/Core_Community_Data2019.csv")
serp <- read.csv("/Users/Marina/Google Drive/02_McLaughlin_80Sites_Organized/Data-Storage-Files/McL_Abiotic-Data.csv")

rm <- c("Logfia spp./Micropus californicus", "Aegilops triuncialis", "Trifolium bifidum/gracilentum", "Torilis arvensis/nodosa", "Hesperolinon spp.", "Cuscuta californica", "Silene gallica", "Erodium botrys", "Erodium brachycarpum", "Trifolium microdon", "Trifolium microcephalum") # rm commonly mis-ID'd species, species that hadnt always been censused (cuscuta) and goatgrass which is treated at the reserve and occasionally pulled from transects; 

# from susan email 8/9/22: "All Geranium should be dissectum"; "Silene gallica may be a mistake"
# from personal investigations: 
## EROBOT & EROBRA: there are plots where most are one, then they switch for one or two years, then go back (e.g., see site 36: EROBRA throughout but then EROBOT in 2010 and 2013). plus notes that say "should be botrys" but then aren't recorded as such
# Microdon and microcephalum: see sites 40, 59, 67 in 2013 for examples, notes saying ID uncertain/not flowering/could be microdon or microcephalum

mcl <- filter(mcl, Species_Name %in% trait[trait$Annual.Perennial == "Annual",]$Species_Name, !(Species_Name %in% rm)) 
mcl <- merge(serp[,1:2], mcl, by = "Site")

mcl[mcl$Species_Name == "Geranium molle",]$Species_Name <- "Geranium dissectum"

mcl$Species_Name <- recode_factor(mcl$Species_Name, 'Lysimachia arvensis' = "Anagallis arvensis", 'Gastridium phleoides' = "Gastridium ventricosum")

mcl <- ddply(mcl, .(Year, Site, Quadrat, Serpentine, Species_Name), summarize, Cover = sum(Cover))

mcl.RA <- mcl %>%
  filter(Year > 2005) %>%
  dplyr::mutate(sum.abundance = sum(Cover, na.rm = T)) %>%
  dplyr::group_by(Species_Name, sum.abundance) %>%
  dplyr::summarize(sp.RA = sum(Cover, na.rm = T)) %>%
  dplyr::group_by(Species_Name) %>%
  dplyr::summarize(RA = sp.RA/sum.abundance)


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

#### metadata ####
mcl$sumPA <- rowSums(mcl[,3:22])

meta.mcl <- mcl %>%
  dplyr::group_by(Species_Name) %>%
  dplyr::summarize(n.plots = length(ID),
            mean.years = mean(sumPA),
            min.years = min(sumPA),
            max.years = max(sumPA))

mcl <- mcl[,23]

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


df <- data.frame(names = NA, obs = NA)

for(n in names(species.list)) {
  df <- data.frame(rbind(df, c(names = n, obs = nrow(species.list[[n]]))))
}


df$obs <- as.numeric(as.character(df$obs))


#### Run HMM ####
mcl.df <- expand.grid(Species_Name = names(species.list), p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA, n.plots = NA, min.years = NA, max.years = NA, mean.years = NA) # empty df to fill with rates


for(j in names(species.list)){ # for each species
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

  mcl.df[mcl.df$Species_Name == j,]$p0 <- EMresult$param$p0
  mcl.df[mcl.df$Species_Name == j,]$g <- EMresult$param$g
  mcl.df[mcl.df$Species_Name == j,]$c <- EMresult$param$c
  mcl.df[mcl.df$Species_Name == j,]$s <- EMresult$param$s
  mcl.df[mcl.df$Species_Name == j,]$r <- EMresult$param$r
  mcl.df[mcl.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
  mcl.df[mcl.df$Species_Name == j,]$n.plots <- nrow(species.list[[j]])
  mcl.df[mcl.df$Species_Name == j,]$mean.years <- mean(rowSums(species.list[[j]]))
  mcl.df[mcl.df$Species_Name == j,]$min.years <- min(rowSums(species.list[[j]]))
  mcl.df[mcl.df$Species_Name == j,]$max.years <- max(rowSums(species.list[[j]]))
}


#### Explore Output ####
trait$Species_Name <- recode_factor(trait$Species_Name, 'Lysimachia arvensis' = "Anagallis arvensis", 'Gastridium phleoides' = "Gastridium ventricosum")

trait$FunGroup <- paste(trait$Native.Exotic, trait$Grass.Forb.Shrub)

meta.mcl <- merge(meta.mcl, trait[,c(1,23)], by = "Species_Name", all.y = F)

#get codes
code <- read.csv("Data/20211001_Full-Species-List.csv")

meta.mcl <- merge(meta.mcl, code[,c(2,6)], by.x = "Species_Name", by.y = "Species", all.x = T)

meta.mcl$new.code <- ifelse(is.na(meta.mcl$new.code), 
               sapply(str_split(meta.mcl$Species_Name, " "), function(x){
      toupper(paste(substring(x, 1, 3), collapse = "")) }),
  meta.mcl$new.code)

write.csv(meta.mcl, "/Users/Marina/Documents/USDA-PostDoc/Projects/Seed-Traits/Data/Long-Term-Datasets/McLaughlin/HMM-meta-mcl.csv", row.names = F)

mcl.df <- merge(mcl.df, trait[,c(1,23)], by = "Species_Name", all.y = F)

# mcl <- merge(mcl, trait[,c(1,23)], by = "Species_Name", all.y = F)
# mcl$site.id <- paste(mcl$Site, mcl$Quadrat, sep = ".")

mcl.RA <- merge(mcl.RA, trait[,c(1,3)], by = "Species_Name", all.x = T, all.y = F)
mcl.RA <- mcl.RA[,3:2]
saveRDS(mcl.RA, "/Users/Marina/Documents/USDA-PostDoc/Projects/Seed-Traits/Data/Long-Term-Datasets/McLaughlin/McL_RA.RDS")

ggplot(mcl.df[mcl.df$iter < 150 & mcl.df$n.plots >= 50,], aes(x = s, y = c, col = FunGroup)) +
  #geom_point() +
  geom_text(aes(label = Code)) + 
  geom_smooth(method = "lm", formula = y ~ x, se = F) +
  theme_bw()

saveRDS(mcl.df, "Scripts/HMMs/McLaughlin/20230704_mcl_HMM_quadrat.RDS")

