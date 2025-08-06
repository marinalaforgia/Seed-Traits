#### Jasper Ridge HMM ####
rm(list=ls())
source("Scripts/HMMs/HMM_Functions.R")

#### Load Libraries ####
library(ggplot2)
library(tidyverse)

#### Prep Data ####
jr <- read.csv("Data/Long-Term-Datasets/Jasper-Ridge/JR_cover_updated.csv")[,-1] #weird first column

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

jr2 <- filter(jr, cover > 0, treatment == "c", !species %in% rm)

# overall mean
jr2$sum.abundance <- sum(jr2$cover)

jr.RA <- jr2 %>%
  #dplyr::mutate(sum.abundance = sum(cover)) %>% # WHY DOESNT THIS WORK
  dplyr::group_by(species, sum.abundance) %>%
  dplyr::summarize(sp.RA = sum(cover)) %>%
  dplyr::group_by(species) %>%
  dplyr::summarize(RA = mean(sp.RA)/sum.abundance)

#unique(jr$Site)
#sites <- c("83B", "88", "92", "89") # for counts
#jr <- filter(jr, Site %in% sites) # for counts
#jr <- filter(jr, treatment == "r")
#jr <- filter(jr, species != "BRHO")
#jr <- filter(jr, species != "PLER")
jr.sum <- pivot_wider(jr[,c(1:3,6)], names_from = "year", values_from = "PA")
jr.sum[is.na(jr.sum)] <- 0
jr.sum$sumPA <- rowSums(jr.sum[,3:39])
jr.sum <- filter(jr.sum, sumPA > 0) 

jr.sum <- filter(jr.sum, !species %in% rm)

#control.plots <- c("c101", "c103", "c105", "c108", "c110", "c112", "c113", "c115", "c117", "c120", "c122", "c124", "c201", "c203", "c205", "c208", "c210", "c212", "c213", "c215", "c217", "c220", "c222", "c224", "c301", "c303", "c305", "c308", "c310", "c312", "c313", "c315", "c317", "c320", "c322", "c324")

control.plots <- c("c102", "c104", "c106", "c107", "c109", "c111", "c114", "c116", "c118", "c119", "c121", "c123", "c202", "c204", "c206", "c207", "c209", "c211", "c214", "c216", "c218", "c219", "c221", "c223", "c302", "c304", "c306", "c307", "c309", "c311", "c314", "c316", "c318", "c319", "c321", "c323")

#gopher.plots <- c("g101", "g103", "g105", "g108", "g110", "g112", "g113", "g115", "g117", "g120", "g122", "g124", "g201", "g203", "g205", "g208", "g210", "g212", "g213", "g215", "g217", "g220", "g222", "g224", "g301", "g303", "g305", "g308", "g310", "g312", "g313", "g315", "g317", "g320", "g322", "g324")

#rabbit.plots <- c("r101", "r103", "r105", "r108", "r110", "r112", "r113", "r115", "r117", "r120", "r122", "r124", "r201", "r203", "r205", "r208", "r210", "r212", "r213", "r215", "r217", "r220", "r222", "r224", "r301", "r303", "r305", "r308", "r310", "r312", "r313", "r315", "r317", "r320", "r322", "r324")

#all <- c(control.plots, rabbit.plots, gopher.plots)
jr.sum <- filter(jr.sum, quadID %in% control.plots)


#### Metadata ####
meta.jr <- jr.sum %>%
  dplyr::group_by(species) %>%
  dplyr::summarize(n.plots = length(quadID),
            mean.years = mean(sumPA),
            min.years = min(sumPA),
            max.years = max(sumPA))

jr.sum <- jr.sum[,-40]

#jr.sum <- filter(jr.sum, QuadratCode %in% gopher.plots)

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
# list; each object is a dataframe of PA data for a single species where each row is a site and each column is a year
jr.df <- expand.grid(Species_Name = names(species.list), p0 = NA, g = NA, c = NA, s = NA, r = NA, iter = NA, n.plots = NA, min.years = NA, max.years = NA, mean.years = NA) # empty df to fill with rates

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
    for (k in 1:5) print(i) # I don't understand what's going on here, why print it 5 times?
    EMresult = EMestimation(X, r = 1) # update the model params, stop when new log likelihood is lower than the old log likelihood plus some precision 
    # fill in results with new params from the best log likelihood model
    p0Results[i] = EMresult$param$p0 # initial seed bank prob
    gResults[i] = EMresult$param$g # germ
    cResults[i] = EMresult$param$c #col
    sResults[i] = EMresult$param$s #surv
    rResults[i] = EMresult$param$r #prod
    llResults[i] = EMresult$ll # log-likelihood
  }

  # here I assumed that the last estimated parameters, were what the model converged on but I'm not sure that's correct?
  jr.df[jr.df$Species_Name == j,]$p0 <- EMresult$param$p0
  jr.df[jr.df$Species_Name == j,]$g <- EMresult$param$g
  jr.df[jr.df$Species_Name == j,]$c <- EMresult$param$c
  jr.df[jr.df$Species_Name == j,]$s <- EMresult$param$s
  jr.df[jr.df$Species_Name == j,]$r <- EMresult$param$r
  jr.df[jr.df$Species_Name == j,]$iter <- which.max(EMresult$lllist)
  jr.df[jr.df$Species_Name == j,]$n.plots <- nrow(species.list[[j]])
  jr.df[jr.df$Species_Name == j,]$mean.years <- mean(rowSums(species.list[[j]]))
  jr.df[jr.df$Species_Name == j,]$min.years <- min(rowSums(species.list[[j]]))
  jr.df[jr.df$Species_Name == j,]$max.years <- max(rowSums(species.list[[j]]))
}


# remove species that didn't converge in 150 steps
#jr.df <- filter(jr.df, iter < 100)

#### Explore Output ####

trait <- read.csv("Data/Long-Term-Datasets/Jasper-Ridge/JR_Traits.csv")
trait$Code.4 <- toupper(trait$Code.4)
trait$FunGroup <- paste(trait$Status, trait$Habit)

jr.df <- merge(jr.df, unique(trait[,c(1,2,3,23)]), by.x = "Species_Name", by.y = "Code.4", all.x = T, all.y = F)
jr.df <- jr.df[,c(12:14,2:11)]
colnames(jr.df)[1] <- "Species_Name"
saveRDS(jr.df, "Scripts/HMMs/Jasper-Ridge/jr_HMM_05m_control2.RDS")

jr.df2 <- readRDS("Scripts/HMMs/Jasper-Ridge/jr_HMM_05m_control.RDS")

jr.df <- rbind(jr.df, jr.df2)

test <- jr.df %>%
  group_by(Species_Name, Code, FunGroup) %>%
  summarize(across(.cols = everything(), mean))

ggplot(test[test$iter < 150,], aes(x = s, y = c)) +
  #geom_point() +
  geom_text(aes(label = Code)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_classic()

saveRDS(test, "Scripts/HMMs/Jasper-Ridge/jr_HMM_05m_control-avg.RDS")

## merge with RA
jr.RA <- merge(jr.RA, trait[,2:3], by.x = "species", by.y = "Code.4", all.x = T, all.y = F)
jr.RA <- jr.RA[,c(3,2)]

saveRDS(jr.RA, "Data/Long-Term-Datasets/Jasper-Ridge/JR_RA.rds")

## merge metadata

meta.jr <- merge(meta.jr, unique(trait[,c(1,2,3,23)]), by.x = "species", by.y = "Code.4", all.x = T, all.y = F)
meta.jr <- meta.jr[,c(6:8,2:5)]
write.csv(meta.jr, "Data/Long-Term-Datasets/Jasper-Ridge/HMM-meta-JR-A.csv", row.names = F)

meta.jrA <- read.csv("Data/Long-Term-Datasets/Jasper-Ridge/HMM-meta-JR-A.csv")
meta.jrB <- read.csv("Data/Long-Term-Datasets/Jasper-Ridge/HMM-meta-JR-B.csv")

meta.jr <- rbind(meta.jrA, meta.jrB)
colnames(meta.jr)[1] <- "Species_Name"

meta.jr <- meta.jr %>%
  group_by(Species_Name, Code, FunGroup) %>%
  summarize(across(n.plots:max.years, ~mean(.x)))

write.csv(meta.jr, "Data/Long-Term-Datasets/Jasper-Ridge/HMM-meta-JR-avg.csv", row.names = F)
