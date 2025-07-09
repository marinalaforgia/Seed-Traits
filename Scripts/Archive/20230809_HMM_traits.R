#### HMM System Comparison ####
rm(list=ls())

# what I SHOULD do is run the models multiple times, extract se for each species and throw out the ones that jump around, but first i need to get either the cluster or jetstream to work

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

#### Prep Long term Data ####
species <- read.csv("Data/20211001_Full-Species-List.csv")

# Carrizo #
## TWO OPTIONS: remove first year or not
## 20% of plots = 7.2
CP.df <- readRDS("Scripts/HMMs/Carrizo-Plain/CP_HMM_none-control_ungrazed-GKR_2008.RDS")
#CP.df <- filter(CP.df, Code != "TRIGRA") # removing TRIGRA because... it's estimate was crazy ?? 

# Jasper Ridge #
## 20% of plots = 7.2 
# ran two HMMs and took the average
jr.df <- readRDS("Scripts/HMMs/Jasper-Ridge/jr_HMM_05m_control-avg.RDS") 

# Sonoran Desert #
## 20% of plots = min 8.8
sd.df <- readRDS("Scripts/HMMs/Sonoran-Desert/SD_HMM_2001_01.RDS")

# McLaughlin #
## 20% of plots = min 80
mcl.df <- readRDS("Scripts/HMMs/McLaughlin/20221212_mcl_HMM_quadrat.RDS")[,-14]
mcl.df <- filter(mcl.df, n.plots >= 80)

# Portal #
## 20% of plots = min 13
port.df <- readRDS("Scripts/HMMs/Portal/Port_HMM_control_winter_2009.RDS")

# Add sites #
sd.df$site <- "Sonoran"
port.df$site <- "Portal"
mcl.df$site <- "McLaughlin"
CP.df$site <- "Carrizo"
jr.df$site <- "Jasper Ridge"

meta <- read.csv("Data/Long-Term-Datasets/Datasets_metadata.csv")

# bind together #
full <- rbind(CP.df, jr.df, mcl.df, sd.df, port.df)

rm(CP.df, jr.df, mcl.df, sd.df, port.df)

full <- merge(full, meta[,c(1,2,19)], by.x = "site", by.y = "Dataset")

#20% of plots was way too low for every site except mclaughlin, trying a minimum of 15. ideally more but the more stringent the cutoff, the less species we have to work with and the more McLaughlin just dominates
full <- filter(full, iter < 150, n.plots >= 15)

#### Prep Seed Traits ####
traits <- read.csv("Data/20230801_Seed-Traits_clean_site.csv")

species.AI <- read.csv("Data/20230530_Seeds_All-Accessions.csv")

# AMSMEN has a suspiciously low prop.C compared to other Amsinckia species, replace with mean
traits[traits$Species == "Amsinckia menziesii",]$prop.C <- mean(traits$prop.C, na.rm = T)

traits <- filter(traits, new.code != "CUCPAL") # remove because it's so different, doesnt represent the range of annuals well
traits[is.na(traits$prop.C),]$prop.C <- mean(traits$prop.C, na.rm = T) # HERHIR
traits[is.na(traits$prop.N),]$prop.N <- mean(traits$prop.N, na.rm = T) # HERHIR
traits[is.na(traits$coat.perm.perc),]$coat.perm.perc <- mean(traits$coat.perm.perc, na.rm = T) # HERHIR

traits$FunGroup <- paste(traits$nat.inv, traits$group, sep = " ")

traits$FunGroup <- ifelse(traits$FunGroup == "native forb", "Native Forb", ifelse(traits$FunGroup == "native grass", "Native Grass", traits$FunGroup))

traits$FunGroup <- ifelse(traits$FunGroup == "invasive forb", "Exotic Forb", ifelse(traits$FunGroup == "invasive grass", "Exotic Grass", traits$FunGroup))

traits <- traits %>%
  group_by(Species, new.code, family, group, nat.inv, FunGroup, appendage.type, appendage, disp.cat.all, disp.cat.nat)  %>%
  dplyr::summarize(across(morph.mass.mg:carb.prop.kew, ~ mean(.x, na.rm = TRUE)))

traits.NS <- traits

# Normalize all traits
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

# need to go back and take these species out of MCL see susan email about ID confusion
hmm.trait <- filter(hmm.trait, Code != "GERDIS", Code != "GERMOLL", Code != "SILGAL")

hmm.trait.forb <- filter(hmm.trait.forb, Code != "GERDIS", Code != "GERMOLL", Code != "SILGAL")

hmm.trait$site <- factor(hmm.trait$site, levels = c("Sonoran", "Portal", "Carrizo", "Jasper Ridge", "McLaughlin"))

hmm.trait.forb$site <- factor(hmm.trait.forb$site, levels = c("Sonoran", "Portal", "Carrizo", "Jasper Ridge", "McLaughlin"))

#traits.site.f$site <- factor(traits.site.f$site, levels = c("Sonoran", "Portal", "Carrizo", "Jasper Ridge", "McLaughlin"))

#### K-means clustering ####

### MCL ###
# set.seed(123)
# fviz_nbclust(hmm.trait[hmm.trait$site == "McLaughlin", c(7,8)], kmeans, method = "wss") #3

set.seed(123)
(k.mcl <- kmeans(hmm.trait[hmm.trait$site == "McLaughlin",c(7,8)], 3, nstart = 25))
# 1: medium S - medium T
# 2: high S - low T
# 3: low S - high T 
k.mcl <- data.frame(Code = hmm.trait[hmm.trait$site == "McLaughlin",]$Code, k.means = k.mcl$cluster, site = "McLaughlin")

k.mcl$k.cat <- ifelse(k.mcl$k.means == 1, "medium S - medium T",
                     ifelse(k.mcl$k.means == 2, "high S - low T",
                            "low S - high T"))

### Carrizo ###
set.seed(123)
fviz_nbclust(hmm.trait[hmm.trait$site == "Carrizo", c(7,8)], kmeans, method = "wss") #3

set.seed(123)
(k.cp <- kmeans(hmm.trait[hmm.trait$site == "Carrizo",c(7,8)], 3, nstart = 25))
k.cp <- data.frame(Code = hmm.trait[hmm.trait$site == "Carrizo",]$Code, k.means = k.cp$cluster, site = "Carrizo")
# 1: medium S - medium T
# 2: low S - high T
# 3: high S - low T

# without 2007
k.cp$k.cat <- ifelse(k.cp$k.means == 1, "low S - high T",
                     ifelse(k.cp$k.means == 2, "high S - low T",
                            "medium S - medium T"))
# with 2007
# k.cp$k.cat <- ifelse(k.cp$k.means == 1, "high S - low T",
#                      ifelse(k.cp$k.means == 2, "medium S - medium T",
#                             "low S - high T"))

### Portal ###
set.seed(123)
fviz_nbclust(hmm.trait[hmm.trait$site == "Portal",c(7,8)], kmeans, method = "wss") # 3

set.seed(123)
(k.portal <- kmeans(hmm.trait[hmm.trait$site == "Portal",c(7,8)], 3, nstart = 25))
# 1: high S - medium T
# 2: low S - high T
# 3: high S - low T

k.portal <- data.frame(Code = hmm.trait[hmm.trait$site == "Portal",]$Code, k.means = k.portal$cluster, site = "Portal")

k.portal$k.cat <- ifelse(k.portal$k.means == 1, "low S - high T",
                         ifelse(k.portal$k.means == 2, "high S - low T", 
                                "medium S - medium T"))

### Sonoran ###

set.seed(123)
fviz_nbclust(hmm.trait[hmm.trait$site == "Sonoran",c(7,8)], kmeans, method = "wss") # 3

set.seed(123)
(k.sd <- kmeans(hmm.trait[hmm.trait$site == "Sonoran",c(7,8)], 3, nstart = 25))

k.sd <- data.frame(Code = hmm.trait[hmm.trait$site == "Sonoran",]$Code, k.means = k.sd$cluster, site = "Sonoran")
# 1: medium S - high T
# 2: high S - low T
# 3: low S - high T
k.sd$k.cat <- ifelse(k.sd$k.means == 1, "high S - low T",
                     ifelse(k.sd$k.means == 2, "medium S - medium T",
                            "low S - high T"))

### Jasper ###
set.seed(123)
fviz_nbclust(hmm.trait[hmm.trait$site == "Jasper Ridge",c(7,8)], kmeans, method = "wss") #3

set.seed(123)
(k.jr <- kmeans(hmm.trait[hmm.trait$site == "Jasper Ridge",c(7,8)], 3, nstart = 25))
k.jr <- data.frame(Code = hmm.trait[hmm.trait$site == "Jasper Ridge",]$Code, k.means = k.jr$cluster, site = "Jasper Ridge")
# 1: medium S - medium T
# 2: high S - low T
# 3: low S - high T

k.jr$k.cat <- ifelse(k.jr$k.means == 1, "high S - low T",
                     ifelse(k.jr$k.means == 2, "low S - high T",
                            "medium S - medium T"))

## bind together ##
k.means <- rbind(k.mcl, k.cp, k.jr, k.portal, k.sd)

hmm.trait <- merge(hmm.trait, k.means, by = c("Code", "site"))
hmm.trait$k.cat <- factor(hmm.trait$k.cat, levels = c("high S - low T", "medium S - medium T", "low S - high T"))

rm(k.cp, k.jr, k.mcl, k.portal, k.sd, k.means)

#### Fig 1 a-e: Trade-off ####
###### Graph ####

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
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
  scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9)) +
  labs(title = "Sonoran Desert", x = "temporal dispersal", y = "spatial dispersal")

fig1b <- ggplot(hmm.trait[hmm.trait$site == "Portal",], aes(x = s, y = c)) +
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
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
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
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
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
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
  labs(title = "Jasper Ridge", x = "temporal dispersal", y = "spatial dispersal")

fig1e <- ggplot(hmm.trait[hmm.trait$site == "McLaughlin",], aes(x = s, y = c)) +
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
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
  labs(title = "McLaughlin", x = "temporal dispersal", y = "spatial dispersal") 

##### Models ####
cor.test(hmm.trait[hmm.trait$site == "Portal",]$s,
         hmm.trait[hmm.trait$site == "Portal",]$c) # marg (0.051)

cor.test(hmm.trait[hmm.trait$site == "McLaughlin",]$s,
         hmm.trait[hmm.trait$site == "McLaughlin",]$c) #sig

cor.test(hmm.trait[hmm.trait$site == "Jasper Ridge",]$s,
         hmm.trait[hmm.trait$site == "Jasper Ridge",]$c) # sig

cor.test(hmm.trait[hmm.trait$site == "Sonoran",]$s,
         hmm.trait[hmm.trait$site == "Sonoran",]$c) # sig

cor.test(hmm.trait[hmm.trait$site == "Carrizo" & hmm.trait$Code!= "TRIGRA",]$s,
         hmm.trait[hmm.trait$site == "Carrizo" & hmm.trait$Code!= "TRIGRA",]$c) #only sig without TRIGRA

#### Fig 1f: Kmean cats ####

k.cat <- hmm.trait %>%
  dplyr::group_by(site) %>%
  dplyr::mutate(total = length(Species)) %>%
  dplyr::group_by(site, k.cat, total) %>%
  dplyr::summarize(by.cat = length(Species)) %>%
  dplyr::mutate(cat.pro = by.cat/total)

k.cat$k.cat <- factor(k.cat$k.cat, levels = c("low S - high T", "medium S - medium T", "high S - low T"))

fig1f <- ggplot(k.cat, aes(x = site, y = cat.pro, group = k.cat, fill = k.cat)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    #axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.line = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  scale_fill_manual(values = c("#7570B3","#D95F02","#1B9E77")) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  scale_x_discrete(labels = c("SD", "Portal", "CP", "JR", "McL")) +
  labs(y = "proportion of community")


fig1 <- ggarrange(fig1a, fig1b, fig1c, fig1d, fig1e, fig1f, ncol = 3, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), align = "hv", common.legend = TRUE, legend = "bottom") +
  bgcolor("white") +
  border(color = "white")

#ggsave("Manuscript/Aridity/Fig1.png", fig1, dpi = 600, width = 9, height = 6)

rm(fig1a, fig1b, fig1c, fig1d, fig1e, fig1f, fig1)

#### PCA ####

##### All ####
trait.col <- c("Species", "new.code", "family", "group", "nat.inv", "FunGroup", "morph.mass.mg", "chem.mass.mg",  "shape", "size.mm", "appendage.type", "appendage", "prop.C", "prop.N", "set.time.mpsec", "height.cm", "coat.perm.perc", "both.thick", "mucilage", "fleshy.end", "starchy.end", "coat.thick.per.width", "wing.loading", "ldd.all", "ldd.natural")

trait.pca.col <- c("chem.mass.mg",  "shape", 
"prop.C", "prop.N", "set.time.mpsec", "size.mm", "ldd.natural", "coat.thick.per.width", "height.cm", "wing.loading", "coat.perm.perc")

trait.pca <- traits[, trait.col]

pca <- prcomp(trait.pca[, trait.pca.col], scale = T)
summary(pca)

trait.pca <- cbind(trait.pca, pca$x[,1:4])

pca2 <- autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "FunGroup") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "All Species")

hmm.trait2 <- merge(hmm.trait[,c(1,2,7,8,51)], trait.pca, by.x = "Code", by.y = "new.code")

##### PCA with kcat ####
tmp <- merge(hmm.trait[,c(1,2,7,8,51)], trait.pca[,-c(26:29)], by.x = "Code", by.y = "new.code", all = T)

pca <- prcomp(tmp[, trait.pca.col], scale = T)
summary(pca)

pca$x[,2] <- -pca$x[,2]
pca$rotation[,2] <- -pca$rotation[,2]

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

test.pca <- ggplot(tmp, aes(x = PC1, y = PC2)) +
  geom_point(color = "lightgrey", size = 1.5) +
  geom_point(data = tmp[!is.na(tmp$k.cat),], aes(x = PC1, y = PC2, colour = k.cat), inherit.aes = F, size = 2) +
  theme_classic() + 
  theme(
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_rect(fill = NA),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  geom_segment(data = pca_load, 
               aes(x = 0, y = 0, 
                   xend = PC1*5,
                   yend = PC2*5),
               arrow = arrow(length = unit(1/2, 'picas'))) +
  annotate('text', x = (pca_load$PC1*6), y = (pca_load$PC2*6),
           label = pca_load$variable,
           size = 4) +
  scale_color_manual(values = c("#7570B3","#D95F02","#1B9E77"), na.value = "lightgrey")


ggsave("Manuscript/Aridity/PCA-test.png", test.pca, dpi = 600, width = 7, height = 5, units = "in")

##### FIG 3a ####
tmp <- merge(full[,c(1,3)], trait.pca[,c(1:25)], by.x = "Code", by.y = "new.code", all.x = F, all.y = T)

tmp$site <- factor(tmp$site, levels = c("McLaughlin", "Jasper Ridge", "Carrizo", "Portal", "Sonoran"))

pca <- prcomp(tmp[, trait.pca.col], scale = T)
summary(pca)

pca$x[,2] <- -pca$x[,2]
pca$rotation[,2] <- -pca$rotation[,2]

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

fig3a <- ggplot(tmp[!is.na(tmp$site),], aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = site)) +
  theme_classic() + 
  theme(
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    legend.text = element_text(size = 10),
    legend.position = c(0.86, 0.82),
    legend.background = element_rect(fill = NA),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  geom_polygon(data = pca_hull,
               aes(fill = site,
                   colour = site),
               alpha = 0.05,
               show.legend = FALSE) +
  geom_segment(data = pca_load, 
               aes(x = 0, y = 0, 
                   xend = PC1*5,
                   yend = PC2*5),
               arrow = arrow(length = unit(1/2, 'picas'))) +
  annotate('text', x = (pca_load$PC1*6), y = (pca_load$PC2*6),
           label = pca_load$variable,
           size = 4) +
  scale_color_manual(values = c("#018571", "#80cdc1", "#b4c8a8", "#da8a5c", "#ca562c")) +
  scale_fill_manual(values = c("#018571", "#80cdc1", "#b4c8a8", "#da8a5c", "#ca562c"))


##### Figure 2a? ####
test.pca2 <- ggplot(tmp, aes(x = PC1, y = PC2)) +
  #geom_point(color = "lightgrey", size = 1.5) +
  #geom_point(data = tmp[!is.na(tmp$site),], aes(x = PC1, y = PC2, colour = site), inherit.aes = F, size = 2) +
  theme_classic() + 
  geom_text(data = tmp[!is.na(tmp$site),], aes(label = Code)) +
  theme(
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_rect(fill = NA),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  geom_segment(data = pca_load,
               aes(x = 0, y = 0,
                   xend = PC1*5,
                   yend = PC2*5),
               arrow = arrow(length = unit(1/2, 'picas'))) +
  annotate('text', x = (pca_load$PC1*6), y = (pca_load$PC2*6),
           label = pca_load$variable,
           size = 4) +
  scale_color_manual(values = c("#018571", "#80cdc1", "#b4c8a8", "#da8a5c", "#ca562c"), na.value = "lightgrey")

ggsave("Manuscript/Aridity/PCA-test2.png", test.pca2, dpi = 600, width = 7, height = 5, units = "in")

#### PC graphs #### 

##### All ####
ggplot(hmm.trait2, aes(x = PC1, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

ggplot(hmm.trait2, aes(x = PC2, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

ggplot(hmm.trait2, aes(x = PC1, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

# well that's interesting 
ggplot(hmm.trait2, aes(x = PC2, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_classic() +
    theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~site, scales = "free_y") +
  labs(y = "temporal dispersal", title = "All Species")

# ##### All HMM ####
# ggplot(hmm.trait, aes(x = PC1, y = c)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "spatial dispersal")
# 
# ggplot(hmm.trait, aes(x = PC2, y = c)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "spatial dispersal")
# 
# ggplot(hmm.trait, aes(x = PC1, y = s)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "temporal dispersal")
# 
# # well that's interesting 
# ggplot(hmm.trait, aes(x = PC2, y = s)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_classic() +
#     theme(
#     axis.line = element_blank(),
#     panel.border = element_rect(linewidth = 1, fill = NA),
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   facet_wrap(~site, scales = "free_y") +
#   labs(y = "temporal dispersal", title = "All Species")

#### Trait graphs #### 

# ##### Spatial Dispersal ####
# 
# # McLaughlin
# ggplot(hmm.trait2, aes(x = ldd.natural, y = c)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "spatial dispersal", x = "dispersal potential")
# 
# # McLaughlin, maybe carrizo but in the opposite direction?
# ggplot(hmm.trait2, aes(x = ldd.all, y = c)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "spatial dispersal", x = "dispersal potential")
# 
# ggplot(hmm.trait2, aes(x = height.cm, y = c)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "spatial dispersal", x = "height") # is height more important in SD??
# 
# ggplot(hmm.trait2, aes(x = set.time.mpsec, y = c)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "spatial dispersal", x = "settling time (mpsec)")
# 
# ggplot(hmm.trait2, aes(x = chem.mass.mg, y = c)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "spatial dispersal", x = "mass")
# 
# ggplot(hmm.trait2, aes(x = morph.mass.mg, y = c)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "spatial dispersal", x = "mass")
# 
# ggplot(hmm.trait2, aes(x = size.mm, y = c)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "spatial dispersal", x = "size")
# 
# ggplot(hmm.trait2, aes(x = wing.loading, y = c)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "spatial dispersal", x = "wing loading")
# 
# ggplot(hmm.trait, aes(x = shape, y = c)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "spatial dispersal", x = "shape")

##### Temporal Dispersal ####
ggplot(hmm.trait, aes(x = shape, y = s)) +
  #geom_text(aes(label = Code)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal") # trend but only mcl is sig

# ggplot(hmm.trait, aes(x = size.mm, y = s)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "temporal dispersal")
# 
# ggplot(hmm.trait2, aes(x = shape, y = s)) +
#   #geom_text(aes(label = Code)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "temporal dispersal") # trend but only mcl is sig
# 
# ggplot(hmm.trait2, aes(x = size.mm, y = s)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "temporal dispersal")
# 
# ggplot(hmm.trait, aes(x = prop.N, y = s)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "temporal dispersal")
# 
# ggplot(hmm.trait, aes(x = coat.thick.per.width, y = s)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "temporal dispersal", x = "seed coat thickness (log scale)")
# 
# 
# ggplot(hmm.trait2, aes(x = both.thick, y = s)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "temporal dispersal", x = "seed coat thickness")
# 
# ggplot(hmm.trait, aes(x = prop.C, y = s)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "temporal dispersal")
# 
# ggplot(hmm.trait, aes(x = coat.perm.perc, y = s)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "temporal dispersal")
# 

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
  geom_errorbar(aes(ymax = FRic.mean + FRic.se, ymin = FRic.mean - FRic.se, col = site), width = 0) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
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
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymax = FDis.mean + FDis.se, ymin = FDis.mean - FDis.se), width = 0.001)

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
  geom_point(aes(col = site), size = 2) +
  geom_smooth(method = "lm", col = "black", se = F) +
  geom_errorbar(aes(ymin = RaoQ.mean - RaoQ.se, ymax = RaoQ.mean + RaoQ.se, col = site), width = 0) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
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

traits.RA <- traits.RA[,c(1:27,36:37)]

traits.RA <- merge(traits.RA, hmm.trait2[,c(1,2,30,31)], by.y = c("Code", "site"), by.x = c("new.code", "site")) # ALL; PC1 sig, PC2 not 

# CWM and CWV (see Holshof et al 2013)
#traits.RA <- traits.RA[,c(1:11,29:31,12:28,32)]

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
  scale_y_continuous(breaks = c(0,1,2)) +
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
  scale_y_continuous(breaks = c(0,1,2)) +
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

##### Fig 3 ####
fig3 <- ggarrange(fig3a, ggarrange(fig3b, fig3c, nrow = 1, ncol = 2, labels = c("(b)", "(c)")), ggarrange(fig3d, fig3e, nrow = 1, ncol = 2, labels = c("(d)", "(e)")), labels = c("(a)","","", "", ""), nrow = 3, ncol = 1, heights = c(2.5,1,1))

ggsave("Manuscript/Aridity/Fig-3.png", fig3, height = 10, width = 6, dpi = 600, units = "in")
