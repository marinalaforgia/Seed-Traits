# HMM Aridity Gradient
rm(list=ls())

calcSE<-function(x){
  x2<-na.omit(x)
  sd(x2)/sqrt(length(x2))
}

#### Load Libraries ####
library(ggfortify)
library(ggpubr)
library(Rmisc)
library(tidyverse)
library(gridExtra)
library(sjstats)
library(FD)
library(factoextra)
library(RColorBrewer)
library(emmeans)

#### Prep Longterm Data ####
species <- read.csv("Data/20211001_Full-Species-List.csv")

##### CP ####
## TWO OPTIONS: remove first year or not
## max plots = 36
## years = 14 (not including 2007)
CP.df <- readRDS("Scripts/HMMs/Carrizo-Plain/CP-species-boot.RDS")
CP.meta <- read.csv("Data/Long-Term-Datasets/Carrizo-Plain/HMM-meta-CP.csv")

CP.df <- as.data.frame(CP.df)
CP.df <- CP.df %>% 
  dplyr::mutate(across(p0.boot:perc.con, ~as.numeric(.x)))

CP.df <- CP.df[CP.df$s.boot.sd < 0.1 & CP.df$c.boot.sd < 0.1 & CP.df$perc.con > 0.5,]
colnames(CP.df)[1] <- "Code"

CP.df <- CP.df[CP.df$Code != "POASEC",] # remove perennial
CP.df <- CP.df[CP.df$Code != "MUIMAR",] # remove perennial

CP.df <- merge(CP.df[,1:5], CP.meta, by.x = "Code", by.y = "SpeciesCode", all.x = T)
CP.df <- CP.df[,c(6,1,7,2:5,8:11)]

##### JR ####
## max plots = 72 but really two sets of 36
## years = 37
# ran two HMMs and took the average
jr.df <- readRDS("Scripts/HMMs/Jasper-Ridge/JR-species-boot-avg.RDS") 
jr.meta <- read.csv("Data/Long-Term-Datasets/Jasper-Ridge/HMM-meta-JR-avg.csv")

jr.df <- as.data.frame(jr.df)

jr.df <- jr.df %>% 
  dplyr::mutate(across(p0.boot:r.boot, ~as.numeric(.x)))

#jr.df <- jr.df[jr.df$s.boot.sd < 0.1 & jr.df$c.boot.sd < 0.1 & jr.df$perc.con > 0.5,] # prefiltered when boots were averaged, but where did that happen?

jr.df <- merge(jr.df[,1:5], jr.meta, by = "Code", all.x = T)

jr.df <- jr.df[,c(6,1,7,2:5,8:11)]

##### SD ####
## max plots = 44
## years = 16
sd.df <- readRDS("Scripts/HMMs/Sonoran-Desert/SD-species-boot.RDS")
sd.meta <- read.csv("Data/Long-Term-Datasets/Sonoran-Desert/HMM-meta-SD.csv")

sd.df <- as.data.frame(sd.df)

sd.df <- sd.df %>% 
  dplyr::mutate(across(p0.boot:perc.con, ~as.numeric(.x)))

sd.df <- sd.df[sd.df$s.boot.sd < 0.1 & sd.df$c.boot.sd < 0.1 & sd.df$perc.con > 0.5,]

colnames(sd.df)[1] <- "Species.old"
sd.df <- merge(sd.df[,1:5], sd.meta, by.x = "Species.old", by.y = "species")

sd.df <- sd.df[,c(6:8, 2:5, 9:12)]

##### McL ####
### 19 years
### 400 plots max
mcl.df <- readRDS("Scripts/HMMs/McLaughlin/McL-species-boot.RDS")
mcl.meta <- read.csv("Data/Long-Term-Datasets/McLaughlin/HMM-meta-mcl.csv")

mcl.df <- mcl.df %>% 
  dplyr::mutate(across(p0.boot:perc.con, ~as.numeric(.x)))

mcl.df <- mcl.df[mcl.df$s.boot.sd < 0.1 & mcl.df$c.boot.sd < 0.1 & mcl.df$perc.con > 0.5,]

mcl.df$Species_Name <- recode_factor(mcl.df$Species_Name, 'Lysimachia arvensis' = "Anagallis arvensis", 'Gastridium phleoides' = "Gastridium ventricosum")

mcl.df <- merge(mcl.df[,1:5], mcl.meta, by = "Species_Name")

mcl.df <- mcl.df[,c(1,11,10,2:9)]
names(mcl.df)[2] <- "Code"

##### Port ####
## max plots = 64
## years = 22
boot.w <- readRDS("Scripts/HMMs/Portal/port-species-boot-winter.RDS")
port.df <- as.data.frame(boot.w)

# boot.s <- readRDS("Scripts/HMMs/Portal/port-species-boot-summer.RDS")
# boot.s <- as.data.frame(boot.s)
# 
# boot.sw <- readRDS("Scripts/HMMs/Portal/Port_HMM_control_summer-winter_1988-2009.RDS")
# boot.sw <- as.data.frame(boot.sw)

#port.df <- rbind(as.data.frame(boot.w), as.data.frame(boot.s))

port.meta <- read.csv("Data/Long-Term-Datasets/Portal/HMM-meta-port-winter.csv")
#port.meta.s <- read.csv("Data/Long-Term-Datasets/Portal/HMM-meta-port-summer.csv")
#port.meta.sw <- read.csv("Data/Long-Term-Datasets/Portal/HMM-meta-port-summer-winter.csv")

#port.meta <- rbind(port.meta, port.meta.s, port.meta.sw)

port.df <- port.df %>% 
  dplyr::mutate(across(p0.boot:perc.con, ~as.numeric(.x)))

port.df <- port.df[port.df$s.boot.sd < 0.1 & port.df$c.boot.sd < 0.1 & port.df$perc.con > 0.5,]
colnames(port.df)[1] <- "Species.old"
port.df <- merge(port.df[,c(1:5)], port.meta, by.x = "Species.old", by.y = "SpeciesCode")

port.df <- port.df[,c(6:8, 2:5, 9:12)]
#port.df <- filter(port.df, n.plots >= 12) 
#boot.sw <- filter(boot.sw, iter < 150, n.plots > 10)
#boot.sw <- boot.sw[,-c(8,9)]
#colnames(boot.sw)[4:7] <- c("p0.boot", "g.boot", "c.boot", "s.boot")
#port.df <- rbind(port.df, boot.sw)

#port.df$Code <- recode_factor(port.df$Code, 'PANHIRH' = "PANHIR", 'PECPAPP' = "PECPAP", "BOETRII" = "BOETRI")

# Add sites #
sd.df$site <- "Sonoran"
port.df$site <- "Portal"
mcl.df$site <- "McLaughlin"
CP.df$site <- "Carrizo"
jr.df$site <- "Jasper Ridge"

meta <- read.csv("Data/Long-Term-Datasets/Datasets_metadata.csv")

# bind together #
full <- rbind(CP.df, jr.df, mcl.df, sd.df, port.df)

full <- merge(full, meta[,c(1,2,15,23)], by.x = "site", by.y = "Dataset")

#what we need is a combination cut-off of some number of plot-years. I'm going to use mean-years and total number of plots. 

rm(CP.df, jr.df, mcl.df, sd.df, port.df, CP.meta, jr.meta, mcl.meta, sd.meta, port.meta)

colnames(full)[5:8] <- c("p0", "g", "c", "s")

full <- full %>%
  mutate(
    site_extent_factor = case_when(
      study.spatial.extent.km2 < 0.1       ~ 1.0,
      study.spatial.extent.km2 < 1.0       ~ 0.75,
      study.spatial.extent.km2 < 10.0      ~ 0.5,
      TRUE                                 ~ 0.33
    )
  )


# Calculate minimum number of plots-years, toss out any species that occurs in less than 3 plots regardless of how many years is it shows up in
base_threshold <- 50  # Minimum desired plot-years 
full <- full %>%
  mutate(
    improved_min_n.plots = round(base_threshold/(mean.years * site_extent_factor)),
    improved_min_n.plots = ifelse(improved_min_n.plots < 3, 3, improved_min_n.plots),
    include_improved = n.plots >= improved_min_n.plots
  )

full <- full[full$include_improved == T,]

#### Prep Traits ####
traits <- read.csv("Data/20230801_Seed-Traits_clean_site.csv")

species.AI <- read.csv("Data/20230530_Seeds_All-Accessions.csv")

traits <- filter(traits, new.code != "CUCPAL") # remove because it's so different, doesnt represent the range of annuals well
traits[is.na(traits$prop.C),]$prop.C <- mean(traits$prop.C, na.rm = T) # HERHIR ** CP traits didnt merge correctly for propN/C, double check
traits[is.na(traits$prop.N),]$prop.N <- mean(traits$prop.N, na.rm = T) # HERHIR
traits[is.na(traits$coat.perm.perc),]$coat.perm.perc <- mean(traits$coat.perm.perc, na.rm = T) # HERHIR

traits$FunGroup <- paste(traits$nat.inv, traits$group, sep = " ")

traits$FunGroup <- ifelse(traits$FunGroup == "native forb", "Native Forb", ifelse(traits$FunGroup == "native grass", "Native Grass", traits$FunGroup))

traits$FunGroup <- ifelse(traits$FunGroup == "invasive forb", "Exotic Forb", ifelse(traits$FunGroup == "invasive grass", "Exotic Grass", traits$FunGroup))

traits <- traits %>%
  group_by(Species, new.code, family, group, nat.inv, FunGroup, appendage.type, appendage, disp.cat.all, disp.cat.nat)  %>%
  dplyr::summarize(across(morph.mass.mg:carb.prop.kew, ~ mean(.x, na.rm = TRUE)))

traits$new.mass <- ifelse(traits$group == "grass", traits$morph.mass.mg, traits$chem.mass.mg)

traits.NS <- traits

# Normalize all traits

#hist(log(seed.traits$new.mass))
traits$new.mass <- log(traits$new.mass)

#hist(log(traits$wing.loading))
traits$wing.loading <- log(traits$wing.loading)

#hist(sqrt(traits$prop.C))
traits$prop.C <- sqrt(traits$prop.C)

#hist(log(traits$prop.N))
traits$prop.N <- log(traits$prop.N)

#hist(log(traits$coat.perm.perc))
traits$coat.perm.perc <- log(traits$coat.perm.perc)

#hist(log(traits$morph.mass.mg))
traits$morph.mass.mg <- log(traits$morph.mass.mg)

#hist(log(traits$chem.mass.mg))
traits$chem.mass.mg <- log(traits$chem.mass.mg)

#hist(log(traits$size.mm))
traits$size.mm <- log(traits$size.mm)

#hist(traits$set.time.mpsec)

#hist(traits$height.cm)
traits$height.cm <- log(traits$height.cm)

#hist(sqrt(traits$shape))
#traits$shape <- sqrt(traits$shape)

#hist(log(traits$E.S))
traits$E.S <- log(traits$E.S)

#hist(traits$coat.thick.per.width)
traits$coat.thick.per.width <- log(traits$coat.thick.per.width)

#hist(traits$ldd.all)

#hist(traits$ldd.natural)

# traits$AI.group <- ifelse(traits$AI.mean <= 0.2 & traits$AI.75 <= 0.2, "arid", 
#                           ifelse(traits$AI.mean > 0.2 & traits$AI.25 > 0.2, "semi-arid",
#                                  ifelse(traits$AI.mean > 0.2 & traits$AI.25 <= 0.2, "semi-arid preferred", "arid preferred")))
# 
# traits$AI.group <- factor(traits$AI.group, levels = c("arid", "arid preferred", "semi-arid preferred", "semi-arid"))

#traits <- traits[,c(1:26, 41, 27:40, 42)]

#trait.col <- c("Species", "new.code", "family", "group", "nat.inv", "FunGroup", "morph.mass.mg", "chem.mass.mg",  "shape", "size.mm", "appendage.type", "appendage", "prop.C", "prop.N", "set.time.mpsec", "height.cm", "coat.perm.perc", "E.S", "both.thick", "mucilage", "fleshy.end", "starchy.end", "coat.thick.per.width", "ldd.all", "ldd.natural", "AI", "AI.mean", "AI.IQR", "AI.25", "AI.75", "AI.group")

#### Merge HMM & traits ####

hmm.trait <- merge(full, traits, by.x = c("Code", "FunGroup"), by.y = c("new.code", "FunGroup"), all.x = T, all.y = F)

full.forb <- filter(full, FunGroup != "Exotic Grass", FunGroup != "Native Grass")

hmm.trait.forb <- merge(full.forb, traits, by.x = c("Code", "FunGroup"), by.y = c("new.code", "FunGroup"), all.x = T, all.y = F)

hmm.trait$site <- factor(hmm.trait$site, levels = c("Sonoran", "Portal", "Carrizo", "Jasper Ridge", "McLaughlin"))

hmm.trait.forb$site <- factor(hmm.trait.forb$site, levels = c("Sonoran", "Portal", "Carrizo", "Jasper Ridge", "McLaughlin"))

#hmm.trait <- filter(hmm.trait, nat.inv == "native")

#### Temporal vs spatial dispersal across sites ####
# could just say the relative importance of temporal dispersal will increase with aridity

ggplot(hmm.trait.forb, aes(x = site, y = s)) +
  geom_boxplot()

summary(lm(s ~ AI.site, data = hmm.trait.forb))

ggplot(hmm.trait, aes(x = site, y = s)) +
  geom_boxplot()

fig.new <- ggplot(hmm.trait, aes(x = AI.site, y = s, col = site)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  labs(
    x = "Aridity Index",
    y = "Temporal Dispersal"
  ) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  scale_color_manual(values = c(
    "#ca562c", "#da8a5c", "#b4c8a8", "#80cdc1", "#018571"
  ))

summary(lm(s ~ AI.site, data = hmm.trait))

ggsave("Manuscript/Aridity/Fignew.png", fig.new, dpi = 600, width = 5, height = 4)

# ggplot(hmm.trait, aes(x = site, y = c)) +
#   geom_boxplot()
# 
# ggplot(hmm.trait, aes(x = patch.area.m2, y = c, col = AI.site)) +
#   geom_point()

#### K-means ####
### ALL ###
#fviz_nbclust(hmm.trait[, c(7,8)], kmeans, method = "wss")
set.seed(123)
(k.all <- kmeans(hmm.trait[,8], 3, nstart = 25))
k.all <- data.frame(Code = hmm.trait$Code, k.means = k.all$cluster, site = hmm.trait$site, hmm.trait$s, hmm.trait$c)

k.all$k.cat <- ifelse(k.all$k.means == 1, "high temporal dispersal",
                     ifelse(k.all$k.means == 3, "low temporal dispersal",
                            "medium temporal dispersal"))

ggplot(k.all, aes(x = hmm.trait.s, y = hmm.trait.c, col = k.cat)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(aes(col = k.cat, shape = k.cat), size = 2) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    #axis.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
  ) +
  scale_color_manual(values = c("#AA4499","#44AA99","#999933")) +
  facet_wrap(~site, scales = "free_y")


# ### MCL ###
# # set.seed(123)
# # fviz_nbclust(hmm.trait[hmm.trait$site == "McLaughlin", c(7,8)], kmeans, method = "wss") #3
# 
# set.seed(123)
# (k.mcl <- kmeans(hmm.trait[hmm.trait$site == "McLaughlin",c(7,8)], 3, nstart = 25))
# # 1: medium S - medium T
# # 2: high S - low T
# # 3: low S - high T 
# k.mcl <- data.frame(Code = hmm.trait[hmm.trait$site == "McLaughlin",]$Code, k.means = k.mcl$cluster, site = "McLaughlin")
# 
# k.mcl$k.cat <- ifelse(k.mcl$k.means == 1, "high S - low T",
#                      ifelse(k.mcl$k.means == 2, "low S - high T",
#                             "medium S - medium T"))
# 
# ### Carrizo ###
# set.seed(123)
# fviz_nbclust(hmm.trait[hmm.trait$site == "Carrizo", c(7,8)], kmeans, method = "wss") #3
# 
# set.seed(123)
# (k.cp <- kmeans(hmm.trait[hmm.trait$site == "Carrizo",c(7,8)], 3, nstart = 25))
# k.cp <- data.frame(Code = hmm.trait[hmm.trait$site == "Carrizo",]$Code, k.means = k.cp$cluster, site = "Carrizo")
# # 1: medium S - medium T
# # 2: low S - high T
# # 3: high S - low T
# 
# # without 2007
# k.cp$k.cat <- ifelse(k.cp$k.means == 1, "medium S - medium T",
#                      ifelse(k.cp$k.means == 2, "high S - low T",
#                             "low S - high T"))
# # with 2007
# # k.cp$k.cat <- ifelse(k.cp$k.means == 1, "high S - low T",
# #                      ifelse(k.cp$k.means == 2, "medium S - medium T",
# #                             "low S - high T"))
# 
# ### Portal ###
# set.seed(123)
# fviz_nbclust(hmm.trait[hmm.trait$site == "Portal",c(7,8)], kmeans, method = "wss", k.max = 5) # 3
# 
# set.seed(555)
# (k.portal <- kmeans(hmm.trait[hmm.trait$site == "Portal",c(7,8)], 3, nstart = 25))
# # 1: high S - medium T
# # 2: low S - high T
# # 3: high S - low T
# 
# k.portal <- data.frame(Code = hmm.trait[hmm.trait$site == "Portal",]$Code, k.means = k.portal$cluster, site = "Portal")
# 
# k.portal$k.cat <- ifelse(k.portal$k.means == 1, "medium S - medium T",
#                          ifelse(k.portal$k.means == 2, "low S - high T",
#                                 "high S - low T"))
# 
# # k.portal$k.cat <- ifelse(k.portal$k.means == 1, "low S - high T", "medium S - medium T")
# 
# ### Sonoran ###
# 
# set.seed(123)
# fviz_nbclust(hmm.trait[hmm.trait$site == "Sonoran",c(7,8)], kmeans, method = "wss") # 3
# 
# set.seed(123)
# (k.sd <- kmeans(hmm.trait[hmm.trait$site == "Sonoran",c(7,8)], 3, nstart = 25))
# 
# k.sd <- data.frame(Code = hmm.trait[hmm.trait$site == "Sonoran",]$Code, k.means = k.sd$cluster, site = "Sonoran")
# # 1: medium S - high T
# # 2: high S - low T
# # 3: low S - high T
# k.sd$k.cat <- ifelse(k.sd$k.means == 1, "low S - high T",
#                      ifelse(k.sd$k.means == 2, "medium S - medium T",
#                             "high S - low T"))
# 
# ### Jasper ###
# set.seed(123)
# fviz_nbclust(hmm.trait[hmm.trait$site == "Jasper Ridge",c(7,8)], kmeans, method = "wss") #3
# 
# set.seed(123)
# (k.jr <- kmeans(hmm.trait[hmm.trait$site == "Jasper Ridge",c(7,8)], 3, nstart = 25))
# k.jr <- data.frame(Code = hmm.trait[hmm.trait$site == "Jasper Ridge",]$Code, k.means = k.jr$cluster, site = "Jasper Ridge")
# # 1: medium S - medium T
# # 2: high S - low T
# # 3: low S - high T
# 
# k.jr$k.cat <- ifelse(k.jr$k.means == 1, "medium S - medium T",
#                      ifelse(k.jr$k.means == 2, "low S - high T",
#                             "high S - low T"))
# 
# ## bind together ##
# k.means <- rbind(k.mcl, k.cp, k.jr, k.portal, k.sd)
# 
# k.cat.all <-merge(hmm.trait, k.all, by = c("Code"))
hmm.trait <- merge(hmm.trait, k.all[,c(1,3,6)], by = c("Code", "site"))
# hmm.trait$k.cat <- factor(hmm.trait$k.cat, levels = c("high S - low T", "medium S - medium T", "low S - high T"))
# 
# rm(k.cp, k.jr, k.mcl, k.portal, k.sd, k.means)

#### Fig 1a-e: trade-off ####
###### Graph ####

k.plot <- ggplot(k.all, aes(x = hmm.trait.s, y = hmm.trait.c)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(aes(col = k.cat, shape = k.cat), size = 2) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    #axis.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
  ) +
  scale_color_manual(values = c("#AA4499","#44AA99","#999933")) +
  facet_wrap(~site, scales = "free_y")



fig1a <- ggplot(hmm.trait[hmm.trait$site == "Sonoran",], aes(x = s, y = c)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(aes(col = k.cat, shape = k.cat), size = 2) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    #axis.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
  ) +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499")) +
  scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9)) +
  labs(title = "Sonoran Desert", x = "temporal dispersal", y = "spatial dispersal")

fig1b <- ggplot(hmm.trait[hmm.trait$site == "Portal",], aes(x = s, y = c)) +
  #geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(aes(col = k.cat, shape = k.cat), size = 2) +
  theme_classic() +
  #geom_text(aes(label = Code)) +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    #axis.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
  ) +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499")) +
  scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9)) +
  labs(title = "Portal", x = "temporal dispersal", y = "spatial dispersal")

fig1c <- ggplot(hmm.trait[hmm.trait$site == "Carrizo",], aes(x = s, y = c)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(aes(col = k.cat, shape = k.cat), size = 2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    #axis.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
  ) +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499")) +
  scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9)) +
  scale_y_continuous(breaks = c(0.3, 0.5, 0.7)) +
  labs(title = "Carrizo Plain", x = "temporal dispersal", y = "spatial dispersal")

fig1d <- ggplot(hmm.trait[hmm.trait$site == "Jasper Ridge",], aes(x = s, y = c)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(aes(col = k.cat, shape = k.cat), size = 2) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    #axis.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
  ) +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9)) +
  scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9)) +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499")) +
  labs(title = "Jasper Ridge", x = "temporal dispersal", y = "spatial dispersal")

fig1e <- ggplot(hmm.trait[hmm.trait$site == "McLaughlin",], aes(x = s, y = c)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  #geom_point(aes(col = k.cat, shape = k.cat), size = 2) +
  geom_text(aes(label = Code)) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    #axis.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
  ) +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499")) +
  labs(title = "McLaughlin", x = "temporal dispersal", y = "spatial dispersal") 

##### Models ####
cor.test(hmm.trait[hmm.trait$site == "Portal",]$s,
         hmm.trait[hmm.trait$site == "Portal",]$c) # NS

cor.test(hmm.trait[hmm.trait$site == "McLaughlin",]$s,
         hmm.trait[hmm.trait$site == "McLaughlin",]$c) #sig

cor.test(hmm.trait[hmm.trait$site == "Jasper Ridge",]$s,
         hmm.trait[hmm.trait$site == "Jasper Ridge",]$c) # sig

cor.test(hmm.trait[hmm.trait$site == "Sonoran",]$s,
         hmm.trait[hmm.trait$site == "Sonoran",]$c) # sig


cor.test(hmm.trait[hmm.trait$site == "Carrizo",]$s,
         hmm.trait[hmm.trait$site == "Carrizo",]$c) 

#### Fig 1f: prop k cat ####

k.cat <- hmm.trait %>%
  dplyr::group_by(site) %>%
  dplyr::mutate(total = length(Code)) %>%
  dplyr::group_by(site, k.cat, total) %>%
  dplyr::summarize(by.cat = length(Code)) %>%
  dplyr::mutate(cat.pro = by.cat/total)

k.cat$k.cat <- factor(k.cat$k.cat, levels = c("low temporal dispersal", "medium temporal dispersal", "high temporal dispersal"))

fig1f <- ggplot(k.cat, aes(x = site, y = cat.pro, group = k.cat, fill = k.cat)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    #axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.line = element_blank(),
    #legend.position = "none",
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth =  1)
    ) +
  scale_fill_manual(values = c("#44AA99","#999933","#AA4499")) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  scale_x_discrete(labels = c("SD", "Portal", "CP", "JR", "McL")) +
  labs(y = "proportion of community")


fig1 <- ggarrange(fig1a, fig1b, fig1c, fig1d, fig1e, fig1f, ncol = 3, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), align = "hv", common.legend = TRUE, legend = "bottom") +
  bgcolor("white") +
  border(color = "white")

ggsave("Manuscript/Aridity/Fig1.png", fig1, dpi = 600, width = 9, height = 6)

rm(fig1a, fig1b, fig1c, fig1d, fig1e, fig1f, fig1)

##### Figure using k.all ####
k.plot.SD <- ggplot(k.all[k.all$site == "Sonoran",], aes(x = hmm.trait.s, y = hmm.trait.c)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(aes(col = k.cat, shape = k.cat), size = 2) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    #axis.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
  ) +
  scale_color_manual(values = c("#AA4499","#44AA99","#999933")) +
  xlim(0,1) +
  labs(title = "Sonoran Desert", x = "temporal dispersal", y = "spatial dispersal")

k.plot.port <- ggplot(k.all[k.all$site == "Portal",], aes(x = hmm.trait.s, y = hmm.trait.c)) +
  #geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(aes(col = k.cat, shape = k.cat), size = 2) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    #axis.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
  ) +
  scale_color_manual(values = c("#AA4499","#999933")) +
  xlim(0,1) +
  labs(title = "Portal", x = "temporal dispersal", y = "spatial dispersal")

k.plot.CP <- ggplot(k.all[k.all$site == "Carrizo",], aes(x = hmm.trait.s, y = hmm.trait.c)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(aes(col = k.cat, shape = k.cat), size = 2) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    #axis.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
  ) +
  scale_color_manual(values = c("#AA4499","#44AA99","#999933")) +
  xlim(0,1) +
  labs(title = "Carrizo", x = "temporal dispersal", y = "spatial dispersal")

k.plot.JR <- ggplot(k.all[k.all$site == "Jasper Ridge",], aes(x = hmm.trait.s, y = hmm.trait.c)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(aes(col = k.cat, shape = k.cat), size = 2) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    #axis.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
  ) +
  scale_color_manual(values = c("#AA4499","#44AA99","#999933")) +
  xlim(0,1) +
  labs(title = "Jasper Ridge", x = "temporal dispersal", y = "spatial dispersal")

k.plot.Mcl <- ggplot(k.all[k.all$site == "McLaughlin",], aes(x = hmm.trait.s, y = hmm.trait.c)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(aes(col = k.cat, shape = k.cat), size = 2) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    #axis.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
  ) +
  scale_color_manual(values = c("#AA4499","#44AA99","#999933")) +
  xlim(0,1) +
  labs(title = "McLaughlin", x = "temporal dispersal", y = "spatial dispersal")

k.cat <- k.all %>%
  dplyr::group_by(site) %>%
  dplyr::mutate(total = length(Code)) %>%
  dplyr::group_by(site, k.cat, total) %>%
  dplyr::summarize(by.cat = length(Code)) %>%
  dplyr::mutate(cat.pro = by.cat/total)

k.cat$k.cat <- factor(k.cat$k.cat, levels = c("low temporal dispersal", "medium temporal dispersal", "high temporal dispersal"))

fig1f <- ggplot(k.cat, aes(x = site, y = cat.pro, group = k.cat, fill = k.cat)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    #axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.line = element_blank(),
    #legend.position = "none",
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth =  1)
  ) +
  scale_fill_manual(values = c("#44AA99","#999933","#AA4499")) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  scale_x_discrete(labels = c("SD", "Portal", "CP", "JR", "McL")) +
  labs(y = "proportion of community")


fig1 <- ggarrange(k.plot.SD, k.plot.port, k.plot.CP, k.plot.JR, k.plot.Mcl, fig1f, ncol = 3, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), align = "hv", common.legend = TRUE, legend = "bottom") +
  bgcolor("white") +
  border(color = "white")

ggsave("Manuscript/Aridity/Fig1v2.png", fig1, dpi = 600, width = 9, height = 6)

#### Model ####
k.cat.wide <- k.cat[,c(1,2,5)] %>%
  pivot_wider(names_from = k.cat, values_from = cat.pro)

k.cat.wide[k.cat.wide$site == "Portal",3] <- 0

colnames(k.cat.wide)[2] <-"high"
colnames(k.cat.wide)[3] <-"low"
colnames(k.cat.wide)[4] <- "med"

# manova_result <- manova(cbind(high, med, low) ~ site, data = k.cat.wide)
# summary(manova_result, "Wilks")

table_data <- tidyr::pivot_wider(
  k.cat[,c(1:4)],
  names_from = k.cat,
  values_from = by.cat,
  values_fill = 0
)

colnames(table_data)[3] <-"high"
colnames(table_data)[4] <-"low"
colnames(table_data)[5] <- "med"

# Is this more appropriate?
low <- glm(cbind(low, total - low) ~ site, data = table_data, family = binomial)
anova(low) # yes

med <- glm(cbind(med, total - med) ~ site, data = table_data, family = binomial)
anova(med) # yes

high <- glm(cbind(high, total - high) ~ site, data = table_data, family = binomial)
anova(high) # yes

summary(emmeans(low, pairwise ~ site, type = "response")) # nope
summary(emmeans(med, pairwise ~ site, type = "response")) # nope
summary(emmeans(high, pairwise ~ site, type = "response")) # YES
 ## proportion of high dispersers is significantly higher in the two most arid sites, sonoran and portal, than in the 3 semi-arid sites

#### Fig 2a: PCA ####

##### PCA (fungroup) ####
trait.col <- c("Species", "new.code", "family", "group", "nat.inv", "FunGroup", "morph.mass.mg", "chem.mass.mg",  "shape", "size.mm", "appendage.type", "appendage", "prop.C", "prop.N", "set.time.mpsec", "height.cm", "coat.perm.perc", "both.thick", "mucilage", "fleshy.end", "starchy.end", "coat.thick.per.width", "wing.loading", "ldd.all", "ldd.natural")

trait.pca.col <- c("chem.mass.mg",  "shape", 
"set.time.mpsec", "size.mm", "ldd.natural", "wing.loading", "height.cm", "coat.perm.perc", "coat.thick.per.width", "prop.C", "prop.N")

trait.pca <- traits[, trait.col]

pca <- prcomp(trait.pca[, trait.pca.col], scale = T)
summary(pca)

pca$rotation # mass drives PC1, shape drives PC2, prop C PC3, prop N PC4
pca$x[,1] <- -pca$x[,1]
pca$rotation[,1] <- -pca$rotation[,1]

pca$rotation # mass drives PC1, shape drives PC2, prop C PC3, prop N PC4
pca$x[,2] <- -pca$x[,2]
pca$rotation[,2] <- -pca$rotation[,2]

trait.pca <- cbind(trait.pca, pca$x[,1:4])

autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = T, col = "FunGroup") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )


hmm.trait2 <- merge(hmm.trait[,c(1,2,7,8,15)], trait.pca, by.x = "Code", by.y = "new.code")
hmm.trait2 <- merge(k.all[,c(1,3,6)], hmm.trait2, by = c("Code", "site"))

hmm.trait2$k.cat <- recode_factor(hmm.trait2$k.cat, 
                                  `high temporal dispersal` = "high",
                                  `medium temporal dispersal` = "med",
                                  `low temporal dispersal` = "low")
##### PCA (grey) ####
tmp <- merge(hmm.trait[,c(1,2,7,8,48)], trait.pca[,-c(26:29)], by.x = "Code", by.y = "new.code", all.y = T)

pca <- prcomp(tmp[, trait.pca.col], scale = T)
summary(pca)

# pca$x[,2] <- -pca$x[,2]
# pca$rotation[,2] <- -pca$rotation[,2]
# 
# pca$x[,1] <- -pca$x[,1]
# pca$rotation[,1] <- -pca$rotation[,1]

tmp <- cbind(tmp, pca$x[,1:4])

pca_load <- 
  as_tibble(pca$rotation, rownames = 'variable') %>% 
  #we can rename the variables so they look nicer on the figure
  mutate(variable = dplyr::recode(variable,
                  'set.time.mpsec' = 'fall speed',
                  'chem.mass.mg' = 'mass',
                  'coat.thick.per.width' = 'coat thick',
                  'prop.N' = 'N',
                  'prop.C' = 'C',
                  'ldd.natural' = 'dispersal',
                  'size.mm' = 'size',
                  'coat.perm.perc' = 'coat perm',
                  'wing.loading' = 'wing loading',
                  'height.cm' = 'height'
))

pca.gray <- ggplot(tmp, aes(x = PC1, y = PC2)) +
  geom_point(color = "lightgrey", size = 1.5) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    legend.text = element_text(size = 10),
    legend.position = "none",
    legend.background = element_rect(fill = NA),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  lims(x = c(-7,7)) +
  scale_y_continuous(breaks = c(-4,0,4), limits = c(-6.5,7.5)) +
  geom_segment(data = pca_load,
               aes(x = 0, y = 0,
                   xend = PC1*11,
                   yend = PC2*11),
               arrow = arrow(length = unit(1/2, 'picas'))) +
  geom_segment(data = pca_load,
               aes(x = 0, y = 0,
                   xend = -PC1*11,
                   yend = -PC2*11), linetype = 2)

#ggsave("Manuscript/Aridity/pca-grey.png", dpi = 600, width = 6, height = 5, units = "in")

#### Fig 3a: PCA ####
tmp <- merge(full[,c(1,3)], trait.pca[,c(1:25)], by.x = "Code", by.y = "new.code", all.x = F, all.y = T)

tmp$site <- factor(tmp$site, levels = c("McLaughlin", "Jasper Ridge", "Carrizo", "Portal", "Sonoran"))

pca <- prcomp(tmp[, trait.pca.col], scale = T)
summary(pca)

pca$x[,2] <- -pca$x[,2]
pca$rotation[,2] <- -pca$rotation[,2]

pca$x[,1] <- -pca$x[,1]
pca$rotation[,1] <- -pca$rotation[,1]

tmp <- cbind(tmp, pca$x[,1:4])

pca_hull <- 
  tmp %>% 
  filter(!is.na(site)) %>%
  group_by(site) %>% 
  slice(chull(PC1, PC2))

# now, we'll just continue to build on our ggplot object
pca_load <- 
  as_tibble(pca$rotation, rownames = 'variable') %>% 
  #we can rename the variables so they look nicer on the figure
  mutate(variable = dplyr::recode(variable,
                  'set.time.mpsec' = 'fall speed',
                  'chem.mass.mg' = 'mass',
                  'coat.thick.per.width' = 'coat thick',
                  'prop.N' = 'N',
                  'prop.C' = 'C',
                  'ldd.natural' = 'dispersal',
                  'size.mm' = 'size',
                  'E.S' = 'E:S',
                  'coat.perm.perc' = 'coat perm',
                  'wing.loading' = 'wing loading',
                  'height.cm' = 'height'
  
))

fig3a <- ggplot(tmp[!is.na(tmp$site),], aes(x = PC1, y = PC2, fill = site, col = site)) +
  #geom_point(aes(colour = site)) +
  #geom_text(aes(label = Code)) +
  theme_classic() + 
  theme(
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    legend.text = element_text(size = 10),
    legend.position = c(0.15, 0.82),
    legend.background = element_rect(fill = NA),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  geom_polygon(data = pca_hull,
               aes(fill = site,
                   colour = site),
               alpha = 0.05,
               show.legend = T) +
  # geom_segment(data = pca_load, 
  #              aes(x = 0, y = 0, 
  #                  xend = PC1*5,
  #                  yend = PC2*5),
  #              arrow = arrow(length = unit(1/2, 'picas'))) +
  #annotate('text', x = (pca_load$PC1*6), y = (pca_load$PC2*6),
           #label = pca_load$variable,
           #size = 4) +
  scale_color_manual(values = c("#018571", "#80cdc1", "#b4c8a8", "#da8a5c", "#ca562c")) +
  scale_fill_manual(values = c("#018571", "#80cdc1", "#b4c8a8", "#da8a5c", "#ca562c"))

#### Fig 2b-f: PC2 & temporal #### 

library(lmerTest)
summary(lmer(s~PC2 + (1|site), data = hmm.trait2)) # PC2 neg covaries with s

summary(lm(s ~ PC2 + AI.site, data = hmm.trait2)) # sig even without mclaughlin

ggplot(hmm.trait2, aes(x = PC2, y = s, col = site)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  #geom_text(aes(label = Code)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA)
  )

a <- ggplot(hmm.trait2[hmm.trait2$site == "Sonoran",], aes(x = PC2, y = s, col = k.cat, shape = k.cat)) +
  #geom_smooth(method = "lm", se = F, linetype = 2, col = "black") + 
  geom_point(size = 1.5) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),    
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "temporal dispersal", title = "Sonoran Desert") +
  scale_color_manual(values = c("#AA4499", "#999933", "#44AA99"))

b <- ggplot(hmm.trait2[hmm.trait2$site == "Portal",], aes(x = PC2, y = s)) +
  geom_smooth(method = "lm", se = F, linetype = 2, col = "black") + 
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "temporal dispersal", title = "Portal") +
  scale_color_manual(values = c("#AA4499","#999933"))

c <- ggplot(hmm.trait2[hmm.trait2$site == "Carrizo",], aes(x = PC2, y = s)) +
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "temporal dispersal", title = "Carrizo") +
  scale_color_manual(values = c("#AA4499", "#999933", "#44AA99"))

d <- ggplot(hmm.trait2[hmm.trait2$site == "Jasper Ridge",], aes(x = PC2, y = s)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "temporal dispersal", title = "Jasper Ridge") +
  scale_color_manual(values = c("#AA4499", "#999933", "#44AA99"))

e <- ggplot(hmm.trait2[hmm.trait2$site == "McLaughlin",], aes(x = PC2, y = s)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  #geom_label(aes(label = Code)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "temporal dispersal", title = "McLaughlin") +
  scale_color_manual(values = c("#AA4499", "#999933", "#44AA99"))

hmm.trait2.sum <- hmm.trait2 %>%
  #filter(site != "Carrizo") %>%
  group_by(k.cat) %>%
  summarize(s.mean = mean(s),
            s.se = calcSE(s),
            c.mean = mean(c),
            c.se = calcSE(c),
            PC1.mean = mean(PC1),
            PC1.se = calcSE(PC1),
            PC2.mean = mean(PC2),
            PC2.se = calcSE(PC2)
            )

f <- ggplot(hmm.trait2.sum, aes(y = s.mean, x = PC2.mean, col = k.cat, shape = k.cat)) +
  geom_smooth(method = "lm", se = F, linetype = 2, col = "black") + 
  geom_point() +
  # geom_text(y = 0.86, x = -0.1004029, label = "a", col = "black") +
  # geom_text(y = 0.59, x = 0.5295523, label = "b", col = "black") +
  # geom_text(y = 0.32, x = 1.7344799, label = "c", col = "black") +
  geom_errorbar(aes(ymin = s.mean - s.se, ymax = s.mean + s.se), width = 0) +
  geom_errorbarh(aes(xmin = PC2.mean - PC2.se, xmax = PC2.mean + PC2.se), height = 0) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title.y = element_blank()
  ) +
  labs(y = "temporal dispersal", x = "PC2", title = "all sites") +
  #scale_x_continuous(limits = c(-1,3), breaks = c(0, 1, 2)) +
  #scale_y_continuous(limits = c(0.2,0.9), breaks = c(0.2,0.4,0.6,0.8)) +
  scale_color_manual(values = c("#AA4499", "#999933", "#44AA99"))

fig2be.site <- ggarrange(a,b,c,d,e,f, nrow = 2, ncol = 3, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), align = "hv", common.legend = TRUE, legend = "bottom") +
  bgcolor("white") +
  border(color = "white")
  
ggsave("Manuscript/Aridity/s-PC2new.png", fig2be.site, dpi = 600, units = "in", width = 8, height = 5)

##### Models #####
par(mfrow = c(2,2))
m.sd <- lm(s ~ PC2, data = hmm.trait2[hmm.trait2$site == "Sonoran",])
m.sd <- lm(s ~ PC2, data = hmm.trait2[hmm.trait2$site == "Sonoran",])
plot(m.sd)
summary(m.sd) 

m.port <- lm(s ~ PC2, data = hmm.trait2[hmm.trait2$site == "Portal",])
plot(m.port)
summary(m.port) 

m.cp <- lm(s ~ PC2, data = hmm.trait2[hmm.trait2$site == "Carrizo",])
plot(m.cp)
summary(m.cp) # p = 0.547

m.jr <- lm(s ~ PC2, data = hmm.trait2[hmm.trait2$site == "Jasper Ridge",])
plot(m.jr)
summary(m.jr) # p = 0.010

m.mcl <- lm(s ~ PC2, data = hmm.trait2[hmm.trait2$site == "McLaughlin",])
plot(m.mcl)
summary(m.mcl) # p < 0.001

#### Fig 2 ####
##### by site ####
Fig2.site <- ggarrange(test.pca2, fig2be.site, ncol = 1, nrow = 2, heights = c(1.2,1))

ggsave("Manuscript/Aridity/fig-PCA-site.jpg", height = 11, width = 7, units = "in", dpi = 600)

##### by k.cat ####
Fig2.site <- ggarrange(test.pca, fig2be.kcat, ncol = 1, nrow = 2, heights = c(1.2,1))

ggsave("Manuscript/Aridity/fig-PCA-kcat.jpg", height = 11, width = 7, units = "in", dpi = 600)

##### by none ####
Fig2.site <- ggarrange(pca.gray, fig2be.none, ncol = 1, nrow = 2, heights = c(1.2,1))

ggsave("Manuscript/Aridity/fig-PCA-grey.jpg", height = 11, width = 7, units = "in", dpi = 600)

#### Trait Convergence ####
mcl <- readRDS("Data/Long-Term-Datasets/McLaughlin/McL_RA-ts.RDS")
mcl$site <- "McLaughlin"
names(mcl)

jr <- readRDS("Data/Long-Term-Datasets/Jasper-Ridge/JR_RA-ts.rds")
jr$site <- "Jasper Ridge"
colnames(jr)[1] <- "Year"
names(jr)

cp <- readRDS("Data/Long-Term-Datasets/Carrizo-Plain/Carrizo_RA-ts.RDS")
cp$site <- "Carrizo"
names(cp)[1] <- "Year"
colnames(cp)[2] <- "quadID"
colnames(cp)[3] <- "Code"
names(cp)

sd <- readRDS("Data/Long-Term-Datasets/Sonoran-Desert/SD_RA-ts.RDS")
sd$site <- "Sonoran"
names(sd)
colnames(sd)[2] <- "quadID"
names(sd)

portal <- readRDS("Data/Long-Term-Datasets/Portal/Portal_RA-ts.RDS")
portal$site <- "Portal"
colnames(portal)[1] <- "Year"
names(portal)
#portal <- portal[portal$Year != 2006,]

df.RA <- rbind(mcl,jr,cp,sd,portal)

df.RA$Code <- recode_factor(df.RA$Code, 
                            FILCAL = "LOGFIL",
                            LYSARV = "ANAARV",
                            PHLGRA = "MICGRA",
                            PLAOVAI = "PLAINS")

df.RA.mean <- df.RA %>%
  group_by(site, Code) %>%
  summarize(RA = mean(RA, na.rm = T))
  
df.RA.wide <- pivot_wider(df.RA, names_from = "Code", values_from = RA)
df.RA.wide <- as.data.frame(df.RA.wide)
row.names(df.RA.wide) <- paste(df.RA.wide$site, df.RA.wide$Year, df.RA.wide$quadID, sep = ".")
#df.RA.wide <- df.RA.wide[df.RA.wide$Year != 2006,] #portal, 1 species in 2006 but no traits for that species
df.RA.wide <- df.RA.wide[,-c(1:3)]
df.RA.wide <- df.RA.wide[,order(colnames(df.RA.wide))]

traits.FD <- traits[,c(2,12:20,24,25,27)]
traits.FD <- as.data.frame(traits.FD)
row.names(traits.FD) <- traits.FD$new.code
traits.FD <- traits.FD[,-1]
traits.FD <- filter(traits.FD, row.names(traits.FD) %in% colnames(df.RA.wide))
df.RA.wide <- df.RA.wide[, colnames(df.RA.wide) %in% row.names(traits.FD)]
df.RA.wide$sums <- rowSums(df.RA.wide, na.rm = T)
df.RA.wide <- filter(df.RA.wide, sums > 0)

nrow(df.RA.wide[df.RA.wide$sums < 0.50,])/nrow(df.RA.wide) #removes 2.6%
nrow(df.RA.wide[df.RA.wide$sums < 0.75,])/nrow(df.RA.wide) #removes 8.7%

df.RA.wide <- filter(df.RA.wide, sums >= 0.75)
df.RA.wide <- df.RA.wide[,-ncol(df.RA.wide)]
species.sum <- names(which(colSums(df.RA.wide, na.rm = T) == 0))
df.RA.wide <- df.RA.wide[,!colnames(df.RA.wide) %in% species.sum]
traits.FD <- traits.FD[!row.names(traits.FD) %in% species.sum,]
traits.FD <- traits.FD[order(row.names(traits.FD)),]
fd.all.ts <- dbFD(as.matrix(traits.FD), as.matrix(df.RA.wide), w.abun = T, stand.x = T, stand.FRic = T, calc.CWM = T) # WHY WERE SO MANY AXES REMOVED


##### FRic ####
tmp <- as.data.frame(fd.all.ts$FRic)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "FRic"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp$site <- factor(tmp$site, levels = c("McLaughlin", "Jasper Ridge", "Carrizo", "Portal", "Sonoran"))

tmp.sum <- tmp %>% 
  group_by(site, AI.site) %>%
  summarize(FRic.mean = mean(FRic, na.rm = T),
            FRic.se = calcSE(FRic))

fig3b <- ggplot(tmp.sum, aes(x = log(AI.site), y = FRic.mean)) +
  geom_smooth(method = "lm", col = "black", se = F) +
  geom_point(aes(col = site), size = 2) +
  #geom_errorbar(aes(ymax = FRic.mean + FRic.se, ymin = FRic.mean - FRic.se, col = site), width = 0) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    legend.position = "none",
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(y = "FRic", x = "Aridity Index (log)") +
  scale_y_continuous(breaks = c(0.15, 0.25)) +
  scale_color_manual(values = c("#018571", "#80cdc1", "#b4c8a8", "#da8a5c", "#ca562c")) 
 
m1 <- lm(FRic.mean ~ AI.site, data = tmp.sum)
summary(m1)

m1 <- lm(FRic ~ AI.site, data = tmp)
summary(m1) #YESSSS

##### FDis ####
tmp <- as.data.frame(fd.all.ts$FDis)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "FDis"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp.sum <- tmp %>% 
  group_by(site, AI.site) %>%
  summarize(FDis.mean = mean(FDis, na.rm = T),
            FDis.se = calcSE(FDis))

ggplot(tmp.sum, aes(x = AI.site, y = FDis.mean)) +
  geom_point() +
  geom_smooth(method = "lm")
  #geom_errorbar(aes(ymax = FDis.mean + FDis.se, ymin = FDis.mean - FDis.se), width = 0.001)

m1 <- lm(FDis ~ AI.site, data = tmp)
summary(m1) #YESSSS

##### RaoQ ####
tmp <- as.data.frame(fd.all.ts$RaoQ)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "RaoQ"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp$site <- factor(tmp$site, levels = c("McLaughlin", "Jasper Ridge", "Carrizo", "Portal", "Sonoran"))

tmp.sum <- tmp %>% 
  group_by(site, AI.site) %>%
  summarize(RaoQ.mean = mean(RaoQ, na.rm = T),
            RaoQ.se = calcSE(RaoQ))

fig3c <- ggplot(tmp.sum, aes(x = log(AI.site), y = RaoQ.mean)) +
  geom_smooth(method = "lm", col = "black", se = F) +
  geom_point(aes(col = site), size = 2) +
  #geom_errorbar(aes(ymin = RaoQ.mean - RaoQ.se, ymax = RaoQ.mean + RaoQ.se, col = site), width = 0) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    legend.position = "none",
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(x = "Aridity Index (log)", y = "Rao's Q") +  
  scale_y_continuous(breaks = c(5,7)) +
  scale_color_manual(values = c("#018571", "#80cdc1", "#b4c8a8", "#da8a5c", "#ca562c")) 

m1 <- lm(RaoQ ~ AI.site, data = tmp)
summary(m1) #YESSSS

m1 <- lm(RaoQ.mean ~ AI.site, data = tmp.sum)
summary(m1) 

ggplot(tmp, aes(x = log(AI.site), y = RaoQ, group = AI.site)) +
  geom_boxplot() +
  geom_smooth(method = "lm") 

#### CWM traits ####
mcl <- readRDS("Data/Long-Term-Datasets/McLaughlin/McL_RA.RDS")
mcl$site <- "McLaughlin"
names(mcl)

jr <- readRDS("Data/Long-Term-Datasets/Jasper-Ridge/JR_RA.rds")
jr$site <- "Jasper Ridge"

cp <- readRDS("Data/Long-Term-Datasets/Carrizo-Plain/Carrizo_RA.RDS")[,1:2]
cp$site <- "Carrizo"
names(cp)[1] <- "Code"

sd <- readRDS("Data/Long-Term-Datasets/Sonoran-Desert/SD_RA.RDS")
sd$site <- "Sonoran"

portal <- readRDS("Data/Long-Term-Datasets/Portal/Portal_RA.RDS")
portal$site <- "Portal"

df.RA <- rbind(mcl,jr,cp,sd,portal)

df.RA$Code <- recode_factor(df.RA$Code, 
                            FILCAL = "LOGFIL",
                            LYSARV = "ANAARV",
                            PHLGRA = "MICGRA",
                            PLAOVAI = "PLAINS")

traits.RA <- merge(traits.NS, df.RA, by.x = "new.code", by.y = "Code", all = F)

traits.RA <- traits.RA[,c(1:27,36:38)]

traits.RA <- merge(traits.RA, hmm.trait2[,c(1:4,30,31)], by.x = c("new.code", "site"), by.y = c("Code", "site")) # ALL; PC1 sig, PC2 not 

# CWM and CWV (see Holshof et al 2013)
#traits.RA <- traits.RA[,c(1:11,29:31,12:28,32)]

# traits.RA.sp <- traits.RA %>%
#   group_by(site, new.code, s, c) %>%
#   summarize(across(morph.mass.mg:PC2, 
#                    list(cwm = ~ sum(.x * RA, na.rm = T),
#                         cwv = ~ sum(RA * I(.x - sum(.x * RA, na.rm = T))^2, na.rm = T)
#                           )))
# 
# ggplot(traits.RA.sp, aes(x = PC2_cwm, y = s)) +
#   geom_text(aes(label = new.code)) +
#   geom_smooth(method = "lm", col = "black") + 
#   #geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free_y") 


traits.RA <- traits.RA %>%
  group_by(site) %>%
  summarize(across(morph.mass.mg:PC2, 
                   list(cwm = ~ sum(.x * RA, na.rm = T),
                        cwv = ~ sum(RA * I(.x - sum(.x * RA, na.rm = T))^2, na.rm = T)
                          )))


traits.RA <- merge(traits.RA, meta[,c(1,2,19:22)], by.x = "site", by.y = "Dataset")

traits.RA$site <- factor(traits.RA$site, levels = c("McLaughlin", "Jasper Ridge", "Carrizo", "Portal", "Sonoran"))

##### PC1 ####
fig3d <- ggplot(traits.RA, aes(y = PC1_cwm, x = log(AI.site))) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  #geom_text(aes(label = site, col = site)) +
  geom_point(aes(col = site), size = 2) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  scale_color_manual(values = c("#018571", "#80cdc1", "#b4c8a8", "#da8a5c", "#ca562c"))  +
  scale_y_continuous(breaks = c(0,1)) +
  labs(x = "Aridity Index (log)", y = "CWM PC1") 

m1 <- lm(PC1_cwm ~ AI.site, data = traits.RA)
summary(m1)

##### PC2 ####
fig3e <- ggplot(traits.RA, aes(y = PC2_cwm, x = log(AI.site))) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  #geom_text(aes(label = site, col = site)) +
  geom_point(aes(col = site), size = 2) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  scale_y_continuous(breaks = c(0.5,1.5)) +
  scale_color_manual(values = c("#018571", "#80cdc1", "#b4c8a8", "#da8a5c", "#ca562c"))  +
  labs(x = "Aridity Index (log)", y = "CWM PC2") 

m1 <- lm(PC2_cwm ~ AI.site, data = traits.RA)
summary(m1)

# CWV PC2
ggplot(traits.RA, aes(y = PC2_cwv, x = AI.site)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index") 

cor.test(traits.RA$PC2_cwv, traits.RA$AI.site) # p = 0.13

#TRAITS
ggplot(traits.RA, aes(y = coat.thick.per.width_cwm, x = AI.site)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index") 

cor.test(traits.RA$coat.thick.per.width_cwm, traits.RA$AI.site)


a <- ggplot(traits.RA, aes(y = both.thick_cwm, x = AI.site)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index", y = "Seed coat thickness")

cor.test(traits.RA$both.thick_cwm, traits.RA$AI.site)

b <- ggplot(traits.RA, aes(y = shape_cwm, x = log(AI.site))) +
  geom_smooth(method = "lm", se = T) + 
  #geom_point(size = 2) +
  labs(x = "Aridity Index") +
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(y = "CWM Shape")

cor.test(traits.RA$shape_cwm, traits.RA$AI.site)


c <- ggplot(traits.RA, aes(y = size.mm_cwm, x = AI.site)) +
   geom_smooth(method = "lm", se = T) + 
  #geom_point(size = 2) +
  labs(x = "Aridity Index") +
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(y = "CWM Size")

cor.test(traits.RA$size.mm_cwm, traits.RA$AI.site)

  
d <- ggplot(traits.RA, aes(y = ldd.natural_cwm, x = AI.site)) +
   geom_smooth(method = "lm", se = T) + 
  #geom_point(size = 2) +
  labs(x = "Aridity Index") +
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(y = "CWM Dispersal Potential")

cor.test(traits.RA$ldd.natural_cwm, log(traits.RA$AI.site))


ggplot(traits.RA, aes(y = ldd.all_cwm, x = log(AI.site))) +
   geom_smooth(method = "lm", se = T) + 
  #geom_point(size = 2) +
  labs(x = "Aridity Index") +
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(y = "CWM Dispersal Potential")


e <- ggplot(traits.RA, aes(y = prop.N_cwm, x = AI.site)) +
   geom_smooth(method = "lm", se = T) + 
  #geom_point(size = 2) +
  labs(x = "Aridity Index") +
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(y = "CWM Prop N")

cor.test(traits.RA$prop.N_cwm, traits.RA$AI.site)

f <- ggplot(traits.RA, aes(y = prop.C_cwm, x = log(AI.site))) +
   geom_smooth(method = "lm", se = T) + 
  #geom_point(size = 2) +
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index", y = "CWM Prop C")

cor.test(traits.RA$prop.C_cwm, traits.RA$AI.site)

g <- ggplot(traits.RA, aes(y = set.time.mpsec_cwm, x = log(AI.site))) +
   geom_smooth(method = "lm", se = T) + 
  #geom_point(size = 2) +
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index", y = "CWM settling time")

ggplot(traits.RA, aes(y = wing.loading_cwm, x = log(AI.site))) +
   geom_smooth(method = "lm", se = T) + 
  #geom_point(size = 2) +
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index", y = "CWM wing loading")

cor.test(traits.RA$set.time.mpsec_cwm, traits.RA$AI.site)

h <- ggplot(traits.RA, aes(y = chem.mass.mg_cwm, x = AI.site)) +
   geom_smooth(method = "lm", se = T) + 
  #geom_point(size = 2) +
  labs() +
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(y = "CWM mass (mg)", x = "Aridity Index") 

cor.test(traits.RA$chem.mass.mg_cwm, traits.RA$AI.site)

cwm.fig <- ggarrange(a,c,b,d,e,h, ncol = 3, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"))
cwm.fig

#### Fig 3: FULL ####
fig3 <- ggarrange(pca.gray, 
                  ggarrange(fig3a, 
                            ggarrange(fig3b, fig3c, fig3d, fig3e, nrow = 2, ncol = 2, labels = c("(c)", "(d)", "(e)", "(f)"), align = "hv"), labels = c("(b)","","", "", ""), nrow = 2, ncol = 1, heights = c(1.25,1)), labels = c("(a)","","", "", ""), widths = c(1.5, 1))

ggsave("Manuscript/Aridity/Fig-3.png", fig3, height = 7, width = 12, dpi = 600, units = "in")

#### Fig S1: PC1 & spatial ####

a <- ggplot(hmm.trait2[hmm.trait2$site == "Sonoran",], aes(x = PC1, y = c, col = k.cat, shape = k.cat)) +
  #geom_smooth(method = "lm", se = F, linetype = 2, col = "black") + 
  geom_point(size = 1.5) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),    
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  #scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "spatial dispersal", title = "Sonoran Desert") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

b <- ggplot(hmm.trait2[hmm.trait2$site == "Portal",], aes(x = PC1, y = c)) +
  #geom_smooth(method = "lm", se = F, linetype = 2, col = "black") + 
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  #scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "spatial dispersal", title = "Portal") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

c <- ggplot(hmm.trait2[hmm.trait2$site == "Carrizo",], aes(x = PC1, y = c)) +
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  #scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "spatial dispersal", title = "Carrizo") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

d <- ggplot(hmm.trait2[hmm.trait2$site == "Jasper Ridge",], aes(x = PC1, y = c)) +
  #geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  #scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "spatial dispersal", title = "Jasper Ridge") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

e <- ggplot(hmm.trait2[hmm.trait2$site == "McLaughlin",], aes(x = PC1, y = c)) +
  #geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  #scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "spatial dispersal", title = "McLaughlin") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

f <- ggplot(hmm.trait2.sum, aes(y = c.mean, x = PC1.mean, col = k.cat, shape = k.cat)) +
  #geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  #geom_text(y = 0.86, x = -0.1004029, label = "a", col = "black") +
  #geom_text(y = 0.59, x = 0.5295523, label = "b", col = "black") +
  #geom_text(y = 0.32, x = 1.7344799, label = "c", col = "black") +
  geom_errorbar(aes(ymin = c.mean - c.se, ymax = c.mean + c.se), width = 0) +
  geom_errorbarh(aes(xmin = PC1.mean - PC1.se, xmax = PC1.mean + PC1.se), height = 0) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title.y = element_blank()
  ) +
  labs(y = "spatial dispersal", x = "PC1", title = "all sites") +
  #scale_x_continuous(limits = c(-1,3), breaks = c(0, 1, 2)) +
  #scale_y_continuous(limits = c(0.2,0.9), breaks = c(0.2,0.4,0.6,0.8)) +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))


figs1.site <- ggarrange(a,b,c,d,e,f, nrow = 2, ncol = 3, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), align = "hv", common.legend = TRUE, legend = "bottom") +
  bgcolor("white") +
  border(color = "white")
  
ggsave("Manuscript/Aridity/c-PC1.png", figs1.site, dpi = 600, units = "in", width = 8, height = 5)

#### Fig S2: PC1 & temporal ####
a <- ggplot(hmm.trait2[hmm.trait2$site == "Sonoran",], aes(x = PC1, y = s, col = k.cat, shape = k.cat)) +
  #geom_smooth(method = "lm", se = F, linetype = 2, col = "black") + 
  geom_point(size = 1.5) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),    
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  #scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "temporal dispersal", title = "Sonoran Desert") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

b <- ggplot(hmm.trait2[hmm.trait2$site == "Portal",], aes(x = PC1, y = s)) +
  #geom_smooth(method = "lm", se = F, linetype = 2, col = "black") + 
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  #scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "temporal dispersal", title = "Portal") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

c <- ggplot(hmm.trait2[hmm.trait2$site == "Carrizo",], aes(x = PC1, y = s)) +
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  #scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "temporal dispersal", title = "Carrizo") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

d <- ggplot(hmm.trait2[hmm.trait2$site == "Jasper Ridge",], aes(x = PC1, y = s)) +
  #geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  #scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "temporal dispersal", title = "Jasper Ridge") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

e <- ggplot(hmm.trait2[hmm.trait2$site == "McLaughlin",], aes(x = PC1, y = s)) +
  #geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  #scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "temporal dispersal", title = "McLaughlin") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

f <- ggplot(hmm.trait2.sum, aes(y = s.mean, x = PC1.mean, col = k.cat, shape = k.cat)) +
  geom_point() +
  geom_smooth(method = "lm", )
  geom_errorbar(aes(ymin = s.mean - s.se, ymax = s.mean + s.se), width = 0) +
  geom_errorbarh(aes(xmin = PC1.mean - PC1.se, xmax = PC1.mean + PC1.se), height = 0) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title.y = element_blank()
  ) +
  labs(y = "temporal dispersal", x = "PC1", title = "all sites") +
  #scale_x_continuous(limits = c(-1,3), breaks = c(0, 1, 2)) +
  #scale_y_continuous(limits = c(0.2,0.9), breaks = c(0.2,0.4,0.6,0.8)) +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))


figs2.site <- ggarrange(a,b,c,d,e,f, nrow = 2, ncol = 3, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), align = "hv", common.legend = TRUE, legend = "bottom") +
  bgcolor("white") +
  border(color = "white")
  
ggsave("Manuscript/Aridity/s-PC1.png", figs2.site, dpi = 600, units = "in", width = 8, height = 5)

#### Fig S3: PC2 & spatial ####
a <- ggplot(hmm.trait2[hmm.trait2$site == "Sonoran",], aes(x = PC2, y = c, col = k.cat, shape = k.cat)) +
  #geom_smooth(method = "lm", se = F, linetype = 2, col = "black") + 
  geom_point(size = 1.5) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),    
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "spatial dispersal", title = "Sonoran Desert") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

b <- ggplot(hmm.trait2[hmm.trait2$site == "Portal",], aes(x = PC2, y = c)) +
  #geom_smooth(method = "lm", se = F, linetype = 2, col = "black") + 
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "spatial dispersal", title = "Portal") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

c <- ggplot(hmm.trait2[hmm.trait2$site == "Carrizo",], aes(x = PC2, y = c)) +
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "spatial dispersal", title = "Carrizo") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

d <- ggplot(hmm.trait2[hmm.trait2$site == "Jasper Ridge",], aes(x = PC2, y = c)) +
  #geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "spatial dispersal", title = "Jasper Ridge") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

e <- ggplot(hmm.trait2[hmm.trait2$site == "McLaughlin",], aes(x = PC2, y = c)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(size = 1.5, aes(col = k.cat, shape = k.cat)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(-2,6), breaks = c(-2, 0, 2, 4, 6)) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(y = "spatial dispersal", title = "McLaughlin") +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))

f <- ggplot(hmm.trait2.sum, aes(y = c.mean, x = PC2.mean, col = k.cat, shape = k.cat)) +
  #geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  #geom_text(y = 0.86, x = -0.1004029, label = "a", col = "black") +
  #geom_text(y = 0.59, x = 0.5295523, label = "b", col = "black") +
  #geom_text(y = 0.32, x = 1.7344799, label = "c", col = "black") +
  geom_errorbar(aes(ymin = c.mean - c.se, ymax = c.mean + c.se), width = 0) +
  geom_errorbarh(aes(xmin = PC2.mean - PC2.se, xmax = PC2.mean + PC2.se), height = 0) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title.y = element_blank()
  ) +
  labs(y = "spatial dispersal", x = "PC2", title = "all sites") +
  #scale_x_continuous(limits = c(-1,3), breaks = c(0, 1, 2)) +
  #scale_y_continuous(limits = c(0.2,0.9), breaks = c(0.2,0.4,0.6,0.8)) +
  scale_color_manual(values = c("#44AA99", "#999933", "#AA4499"))


figs3.site <- ggarrange(a,b,c,d,e,f, nrow = 2, ncol = 3, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), align = "hv", common.legend = TRUE, legend = "bottom") +
  bgcolor("white") +
  border(color = "white")
  
ggsave("Manuscript/Aridity/c-PC2.png", figs3.site, dpi = 600, units = "in", width = 8, height = 5)


cor.test(hmm.trait2[hmm.trait2$site == "Sonoran",]$shape, hmm.trait2[hmm.trait2$site == "Sonoran",]$s) # 0.15

cor.test(hmm.trait2[hmm.trait2$site == "Sonoran",]$PC2, hmm.trait2[hmm.trait2$site == "Sonoran",]$s) # 0.12

cor.test(hmm.trait2[hmm.trait2$site == "Portal",]$PC2, hmm.trait2[hmm.trait2$site == "Portal",]$s)

cor.test(hmm.trait2[hmm.trait2$site == "Carrizo",]$PC2, hmm.trait2[hmm.trait2$site == "Carrizo",]$s) # 0.56

cor.test(hmm.trait2[hmm.trait2$site == "Jasper Ridge",]$shape, hmm.trait2[hmm.trait2$site == "Jasper Ridge",]$c) # 0.01

cor.test(hmm.trait2[hmm.trait2$site == "McLaughlin",]$PC2, hmm.trait2[hmm.trait2$site == "McLaughlin",]$c) # <0.001

cor.test(hmm.trait2[hmm.trait2$site == "Sonoran" & hmm.trait2$group == "forb",]$size.mm, hmm.trait2[hmm.trait2$site == "Sonoran" & hmm.trait2$group == "forb",]$s) # 0.15

cor.test(hmm.trait2[hmm.trait2$site == "Portal" & hmm.trait2$nat.inv == "native",]$PC2, hmm.trait2[hmm.trait2$site == "Portal" & hmm.trait2$nat.inv == "native",]$s) # 0.02

cor.test(hmm.trait2[hmm.trait2$site == "Portal" & hmm.trait2$group == "forb",]$PC2, hmm.trait2[hmm.trait2$site == "Portal" & hmm.trait2$group == "forb",]$s) # 0.02

cor.test(hmm.trait2[hmm.trait2$site == "Carrizo" & hmm.trait2$group == "forb",]$PC2, hmm.trait2[hmm.trait2$site == "Carrizo" & hmm.trait2$group == "forb",]$s) # 0.03

cor.test(hmm.trait2[hmm.trait2$site == "Jasper Ridge" & hmm.trait2$group == "forb",]$shape, hmm.trait2[hmm.trait2$site == "Jasper Ridge" & hmm.trait2$group == "forb",]$s) # 0.01

cor.test(hmm.trait2[hmm.trait2$site == "McLaughlin" & hmm.trait2$group == "forb",]$shape, hmm.trait2[hmm.trait2$site == "McLaughlin" & hmm.trait2$group == "forb",]$s) # <0.001


#### Random Forest??? ####
# check out Valliere 2020 ecosphere for methods
library(tree)
library(randomForest)
library(rfUtilities)

hmm.trait.NS <- merge(hmm.trait[,c(1,2,7,8)], traits.NS, by.x = "Code", by.y = "new.code")

(s.tree <- tree(s ~ chem.mass.mg + shape + size.mm + coat.perm.perc + coat.thick.per.width + prop.C + prop.N, hmm.trait.NS))

plot(s.tree)

plot(tree(s ~ chem.mass.mg + shape + set.time.mpsec + size.mm + ldd.natural + wing.loading + height.cm + coat.perm.perc + coat.thick.per.width + prop.C + prop.N, hmm.trait.NS))

trait.rf <- randomForest(s ~ chem.mass.mg + shape + size.mm + coat.perm.perc + coat.thick.per.width + prop.C + prop.N, hmm.trait.NS, ntree = 500, importance = T)
print(trait.rf)
summary(trait.rf)
round(importance(trait.rf), 2)

trait.rf <- randomForest(c ~ chem.mass.mg + shape + set.time.mpsec + size.mm + ldd.natural + wing.loading + height.cm, hmm.trait.NS, ntree = 500, importance = T)
print(trait.rf)
summary(trait.rf)
round(importance(trait.rf), 2)

survival <- c("s", "new.mass",  "shape", "size.mm", "prop.C", "prop.N", "coat.perm.perc", "coat.thick.per.width")

colonization <- c("c", "new.mass",  "shape", "set.time.mpsec", "height.cm", "size.mm", "ldd.natural", "wing.loading")
RFMODSEL<-rf.modelSel(x = , y = hmm.trait.NS$s,  imp.scale = "mir", final.model = TRUE)

#### Correlation plots ####


#### Simulation plots ####

##### Mcl ####
mcl.sim.sp <- readRDS("JetStream2/HMM-Sims/McL-sim/McL-sim-species-validation-50.RDS")

mcl.sim.sp <- as.data.frame(mcl.sim.sp)

mcl.sim.sp <- mcl.sim.sp %>% 
  dplyr::mutate(across(sim:iter, ~as.numeric(.x)))

colnames(mcl.sim.sp)[4:9] <- paste(colnames(mcl.sim.sp)[4:9], "est", sep=".")

mcl.sim.sp <- merge(mcl.sim.sp, mcl.df.q, by = "Species_Name")

mcl.rmse <- mcl.sim.sp %>%
  #group_by(sim, patches) %>%
  summarize(s.bias = mean(abs(s - mean(s.est))),
            s.rmse = rmse(s, s.est),
            s.var = var(s, s.est),
            c.bias = mean(abs(c - mean(c.est))),
            c.rmse = rmse(c, c.est),
            c.var = var(c, c.est),
            g.bias = mean(abs(g - mean(g.est))),
            g.rmse = rmse(g, g.est),
            g.var = var(g, g.est)
  )

mcl.sim.sum <- mcl.sim.sp %>%
  dplyr::group_by(Species_Name) %>%
  dplyr::summarize(across(p0.est:s, ~mean(.x)))

a <- ggplot(mcl.sim.sum, aes(y = s.est, x = s)) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_point(size = 0.9) +
  theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
  ) +
  scale_y_continuous(breaks = c(0.3, 0.6)) +
  scale_x_continuous(breaks = c(0.3, 0.6)) +
  labs(x = "real s", y = "mean simulated s")


b <- ggplot(mcl.sim.sum, aes(y = c.est, x = c)) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_point(size = 0.9) +
  theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
  ) +
  labs(x = "real c", y = "mean simulated c")

c <- ggplot(mcl.sim.sum, aes(y = g.est, x = g)) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_point(size = 0.9) +
  theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
  ) +
  labs(x = "real g", y = "mean simulated g")

fig.s2 <- ggarrange(a,b,c, ncol = 3, nrow = 1, labels = c("(a)", "(b)", "(c)"))

ggsave("Manuscript/Figures/FigS2-species-validation.png", fig.s2, height = 2.5, width = 8, units = "in", dpi = 600)