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

trait.pca$hmm <- ifelse(trait.pca$new.code %in% full$Code, "y", "n")

autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "hmm") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(values = c("grey", "red"))


### FIG 3a ####
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


##### Forbs ####
trait.pca.forb <- traits[traits$group == "forb", trait.col]

pca <- prcomp(trait.pca.forb[, trait.pca.col], scale = T)
summary(pca)

trait.pca.forb <- cbind(trait.pca.forb, pca$x[,1:4])

pca.forb2 <- autoplot(pca, x = 1, y = 2, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "FunGroup") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Forbs only")

hmm.trait.forb2 <- merge(hmm.trait[hmm.trait$group == "forb", c(1,2,7,8,51)], trait.pca.forb, by.x = "Code", by.y = "new.code")

fig.pca <- ggarrange(pca.forb2, pca2, widths = c(1,1.3))

#ggsave("Manuscript/Aridity/fig-PCA-all.png", fig.pca, width = 9, height = 4, dpi = 600)

##### All HMM ####
trait.col <- c("Species", "Code", "family", "group", "nat.inv", "FunGroup", "morph.mass.mg", "chem.mass.mg",  "shape", "size.mm", "appendage.type", "appendage", "prop.C", "prop.N", "set.time.mpsec", "height.cm", "coat.perm.perc", "both.thick", "mucilage", "fleshy.end", "starchy.end", "coat.thick.per.width", "wing.loading", "ldd.all", "ldd.natural", "AI", "k.cat")

trait.pca.col <- c("chem.mass.mg",  "shape", 
"prop.C", "prop.N", "set.time.mpsec", "size.mm", "ldd.natural", "coat.thick.per.width", "height.cm", "wing.loading","coat.perm.perc")

trait.pca <- hmm.trait[complete.cases(hmm.trait[,"Species"]), trait.col]
pca <- prcomp(trait.pca[, trait.pca.col], scale = T)
summary(pca)

hmm.trait <- cbind(hmm.trait[complete.cases(hmm.trait[,"Species"]),], pca$x[,1:4])

pca1 <- autoplot(pca, x = 1, y = 2, data = hmm.trait, frame = F, loadings = T, loadings.label = T, label = F, col = "k.cat", size = 3) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
    ) +
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#D95F02"), na.translate = F) +
  stat_ellipse(aes(group = k.cat, col = k.cat)) +
  labs(title = "HMM all species")

autoplot(pca, x = 1, y = 2, data = hmm.trait, frame = F, loadings = T, loadings.label = T, label = F, col = "site", size = 3) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
    ) +
  labs(title = "HMM all species")

m1 <- lm(PC2 ~ site, hmm.trait[hmm.trait$site != "none",])
anova(m1)

m1 <- lm(PC1 ~ site, hmm.trait[hmm.trait$site != "none",])
anova(m1)

m1 <- lm(PC2 ~ k.means, hmm.trait[hmm.trait$site != "none",])
anova(m1)

m1 <- lm(PC1 ~ k.means, hmm.trait[hmm.trait$site != "none",])
anova(m1) # PC1 DOES separate out kmeans, but im sure this is GRASS

pca_hull <- 
  hmm.trait %>% 
  #filter(site != "none") %>%
  group_by(site) %>% 
  slice(chull(PC1, PC2))

# now, we'll just continue to build on our ggplot object
pca_load <- 
  as_tibble(pca$rotation, rownames = 'variable') %>% 
  #we can rename the variables so they look nicer on the figure
  mutate(variable = dplyr::recode(variable,
                  'set.time.mpsec' = 'settling time',
                  'chem.mass.mg' = 'mass',
                  'coat.thick.per.width' = 'coat thickness',
                  'prop.N' = 'N',
                  'prop.C' = 'C',
                  'ldd.natural' = 'dispersal',
                  'size.mm' = 'size',
                  'E.S' = 'E:S',
                  'coat.perm.perc' = 'coat perm',
  
))

ggplot(hmm.trait[hmm.trait$site != "none",], aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = site)) +
  theme_light() + 
  geom_polygon(data = pca_hull,
               aes(fill = site,
                   colour = site),
               alpha = 0.1,
               show.legend = FALSE) +
  geom_segment(data = pca_load, 
               aes(x = 0, y = 0, 
                   xend = PC1*5,
                   yend = PC2*5),
               arrow = arrow(length = unit(1/2, 'picas'))) +
  annotate('text', x = (pca_load$PC1*6), y = (pca_load$PC2*6),
           label = pca_load$variable,
           size = 3.5) 


##### Forbs HMM ####
trait.pca.col <- c("chem.mass.mg",  "shape", "prop.C", "prop.N", "set.time.mpsec", "coat.perm.perc", "size.mm", "ldd.natural", "coat.thick.per.width", "height.cm")

hmm.trait.forb <- merge(hmm.trait.forb, hmm.trait[,c(1,51)], by = "Code", all.x = T, all.y = F)
trait.pca <- hmm.trait.forb[complete.cases(hmm.trait.forb[,"Species"]), trait.col]
pca <- prcomp(trait.pca[, trait.pca.col], scale = T)
summary(pca)

hmm.trait.forb <- cbind(hmm.trait.forb[complete.cases(hmm.trait.forb[,"Species"]),], pca$x[,1:4])

pca.forb1 <- autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "k.cat", size = 3) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
    ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"), na.translate = F) +
  stat_ellipse(aes(group = k.cat, col = k.cat)) +
  labs(title = "HMM forbs")

fig.pca2 <- ggarrange(pca.forb1, pca1, widths = c(1,1.5))

#ggsave("Manuscript/Aridity/fig-PCA-HMM.png", fig.pca2, width = 9, height = 4, dpi = 600)
pca_hull <- 
  hmm.trait.forb %>% 
  group_by(site) %>% 
  slice(chull(PC1, PC2))

# now, we'll just continue to build on our ggplot object
pca_load <- 
  as_tibble(pca$rotation, rownames = 'variable') %>% 
  #we can rename the variables so they look nicer on the figure
  mutate(variable = dplyr::recode(variable,
                  'set.time.mpsec' = 'settling time',
                  'chem.mass.mg' = 'mass',
                  'coat.thick.per.width' = 'coat thickness',
                  'prop.N' = 'N',
                  'prop.C' = 'C',
                  'ldd.natural' = 'dispersal',
                  'size.mm' = 'size',
                  'E.S' = 'E:S',
                  'coat.perm.perc' = 'coat perm',
  
))

ggplot(hmm.trait.forb, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = site)) +
  theme_light() + 
  geom_polygon(data = pca_hull,
               aes(fill = site,
                   colour = site),
               alpha = 0.1,
               show.legend = FALSE) +
  geom_segment(data = pca_load, 
               aes(x = 0, y = 0, 
                   xend = PC1*5,
                   yend = PC2*5),
               arrow = arrow(length = unit(1/2, 'picas'))) +
  annotate('text', x = (pca_load$PC1*6), y = (pca_load$PC2*6),
           label = pca_load$variable,
           size = 3.5) 

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

##### Forbs ####
ggplot(hmm.trait.forb2, aes(x = PC1, y = c)) +
  geom_smooth(method = "lm") +
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

ggplot(hmm.trait.forb2, aes(x = PC2, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

ggplot(hmm.trait.forb2, aes(x = PC1, y = s)) +
  geom_smooth(method = "lm") +
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait.forb2, aes(x = PC2, y = s)) +
  geom_smooth(method = "lm") + 
  #geom_text(aes(label = Code)) + 
  geom_point() +
  theme_classic() +
  geom_vline(xintercept = 0) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA)
  ) +
  facet_wrap(~site, scales = "free_y") +
  labs(y = "temporal dispersal")

##### All HMM ####
ggplot(hmm.trait, aes(x = PC1, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

ggplot(hmm.trait, aes(x = PC2, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

ggplot(hmm.trait, aes(x = PC1, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

# well that's interesting 
ggplot(hmm.trait, aes(x = PC2, y = s)) +
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

##### Forbs HMMM ####
ggplot(hmm.trait.forb, aes(x = PC1, y = c)) +
  geom_smooth(method = "lm") +
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

ggplot(hmm.trait.forb, aes(x = PC1, y = s)) +
  geom_smooth(method = "lm") +
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait.forb, aes(x = PC2, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

ggplot(hmm.trait.forb, aes(x = PC2, y = s)) +
  geom_smooth(method = "lm") + 
  #geom_text(aes(label = Code)) + 
  geom_point() +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_vline(xintercept = 0) +
  facet_wrap(~site, scales = "free_y") +
  labs(y = "temporal dispersal", title = "Forbs Only") # bro are you ok?

#### Trait graphs #### 

##### Spatial Dispersal ####

# McLaughlin
ggplot(hmm.trait2, aes(x = ldd.natural, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "dispersal potential")

# McLaughlin, maybe carrizo but in the opposite direction?
ggplot(hmm.trait2, aes(x = ldd.all, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "dispersal potential")

ggplot(hmm.trait2, aes(x = height.cm, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "height") # is height more important in SD??

ggplot(hmm.trait2, aes(x = set.time.mpsec, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "settling time (mpsec)")

ggplot(hmm.trait2, aes(x = chem.mass.mg, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "mass")

ggplot(hmm.trait2, aes(x = morph.mass.mg, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "mass")

ggplot(hmm.trait2, aes(x = size.mm, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "size")

ggplot(hmm.trait2, aes(x = wing.loading, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "wing loading")

ggplot(hmm.trait, aes(x = shape, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "shape")

##### Temporal Dispersal ####
ggplot(hmm.trait, aes(x = shape, y = s)) +
  #geom_text(aes(label = Code)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal") # trend but only mcl is sig

ggplot(hmm.trait, aes(x = size.mm, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait2, aes(x = shape, y = s)) +
  #geom_text(aes(label = Code)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal") # trend but only mcl is sig

ggplot(hmm.trait2, aes(x = size.mm, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait, aes(x = prop.N, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait, aes(x = coat.thick.per.width, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal", x = "seed coat thickness (log scale)")


ggplot(hmm.trait, aes(x = both.thick, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal", x = "seed coat thickness")

ggplot(hmm.trait, aes(x = prop.C, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait, aes(x = coat.perm.perc, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

##### forbs only #####
ggplot(hmm.trait.forb, aes(x = ldd.natural, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "dispersal potential")

ggplot(hmm.trait.forb, aes(x = set.time.mpsec, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "settling time (mpsec)")

ggplot(hmm.trait.forb, aes(x = size.mm, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait.forb, aes(x = wing.loading, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait.forb, aes(x = shape, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free_y") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait.forb, aes(x = size.mm, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free_y") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait.forb, aes(x = prop.N, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait.forb, aes(x = coat.thick.per.width, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal", x = "seed coat thickness (log scale)")


ggplot(hmm.trait.forb, aes(x = prop.C, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait.forb, aes(x = coat.perm.perc, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait.forb, aes(x = E.S, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

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

#traits.RA <- traits.RA[,c(1:10,36:39,11:27)]

#traits.RA <- merge(traits.RA, hmm.trait[,c(1,2,52,53)], by.y = c("Code", "site"), by.x = c("new.code", "site")) # just HMM; PC1 sig, PC2 not

traits.RA <- merge(traits.RA, hmm.trait2[,c(1,2,30,31)], by.y = c("Code", "site"), by.x = c("new.code", "site")) # ALL; PC1 sig, PC2 not 

#traits.RA <- merge(traits.RA, hmm.trait.forb[,c(1,3,51,52)], by.y = c("Code", "site"), by.x = c("new.code", "site")) #PC1 sig

#traits.RA <- merge(traits.RA, hmm.trait.forb2[,c(1,2,30,31)], by.y = c("Code", "site"), by.x = c("new.code", "site")) # worst


# CWM and CWV (see Holshof et al 2013)
#traits.RA <- traits.RA[,c(1:11,29:31,12:28,32)]

traits.RA <- traits.RA %>%
  group_by(site) %>%
  summarize(across(morph.mass.mg:PC2, 
                   list(cwm = ~ sum(.x * RA, na.rm = T),
                        cwv = ~ sum(RA * I(.x - sum(.x * RA, na.rm = T))^2, na.rm = T)
                          )))


traits.RA <- merge(traits.RA, meta[,c(1,2,19:22)], by.x = "site", by.y = "Dataset")

##### CWM: Aridity index ####

# CWM PC1
pca1 # size and shape
ggplot(traits.RA, aes(y = PC1_cwm, x = log(AI.site))) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index") 

cor.test(traits.RA$PC1_cwm, traits.RA$AI.site) # p = 0.051

ggplot(traits.RA, aes(y = PC1_cwv, x = log(AI.site))) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index") 

cor.test(traits.RA$PC1_cwv, log(traits.RA$AI.site)) # p = 0.13

# CWM PC2
ggplot(traits.RA, aes(y = PC2_cwm, x = AI.site)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index") 

cor.test(traits.RA$PC2_cwm, traits.RA$AI.site) # p = 0.13

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

##### CWM: precip var ####

ggplot(traits.RA, aes(y = coat.thick.per.width_cwm, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")
  
ggplot(traits.RA, aes(y = both.thick_cwm, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")

ggplot(traits.RA, aes(y = shape_cwm, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")

ggplot(traits.RA, aes(y = size.mm_cwm, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")
  
ggplot(traits.RA, aes(y = ldd.natural_cwm, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")

ggplot(traits.RA, aes(y = ldd.all_cwm, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")

ggplot(traits.RA, aes(y = prop.N_cwm, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")

ggplot(traits.RA, aes(y = prop.C_cwm, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")


ggplot(traits.RA, aes(y = set.time.mpsec_cwm, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")

ggplot(traits.RA, aes(y = chem.mass.mg_cwm, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")

##### Predictability ####

ggplot(traits.RA, aes(y = coat.thick.per.width_cwm, x = Colwell.M)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "predictability")
  
ggplot(traits.RA, aes(y = both.thick_cwm, x = Colwell.M)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "predictability")

ggplot(traits.RA, aes(y = shape_cwm, x = Colwell.M)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "predictability")

ggplot(traits.RA, aes(y = size.mm_cwm, x = Colwell.M)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "predictability")
  
ggplot(traits.RA, aes(y = ldd.natural_cwm, x = Colwell.M)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "predictability")

ggplot(traits.RA, aes(y = ldd.all_cwm, x = Colwell.M)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "predictability")

ggplot(traits.RA, aes(y = prop.N_cwm, x = Colwell.M)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "predictability")

ggplot(traits.RA, aes(y = prop.C_cwm, x = Colwell.M)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "predictability")


ggplot(traits.RA, aes(y = set.time.mpsec_cwm, x = Colwell.M)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "predictability")

ggplot(traits.RA, aes(y = chem.mass.mg_cwm, x = Colwell.M)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "predictability")

##### CWV: Aridity index ####
ggplot(traits.RA, aes(y = coat.thick.per.width_cwv, x = AI.site)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index")
  
ggplot(traits.RA, aes(y = both.thick_cwv, x = log(AI.site))) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index")

ggplot(traits.RA, aes(y = shape_cwv, x = log(AI.site))) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index")

ggplot(traits.RA, aes(y = log(size.mm_cwv), x = AI.site)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index")
  
ggplot(traits.RA, aes(y = ldd.natural_cwv, x = log(AI.site))) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index")

ggplot(traits.RA, aes(y = ldd.all_cwv, x = AI.site)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index")

ggplot(traits.RA, aes(y = prop.N_cwv, x = AI.site)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index")

ggplot(traits.RA, aes(y = prop.C_cwv, x = AI.site)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index")


ggplot(traits.RA, aes(y = set.time.mpsec_cwv, x = AI.site)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index")

ggplot(traits.RA, aes(y = log(chem.mass.mg_cwv), x = AI.site)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index")

##### CWV: precip var ####
ggplot(traits.RA, aes(y = coat.thick.per.width_cwv, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")
  
# only one that is significant
ggplot(traits.RA, aes(y = both.thick_cwv, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")

ggplot(traits.RA, aes(y = shape_cwv, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")

ggplot(traits.RA, aes(y = size.mm_cwv, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")
  
ggplot(traits.RA, aes(y = ldd.natural_cwv, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")

ggplot(traits.RA, aes(y = ldd.all_cwv, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")

ggplot(traits.RA, aes(y = prop.N_cwv, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")

ggplot(traits.RA, aes(y = prop.C_cwv, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")


ggplot(traits.RA, aes(y = set.time.mpsec_cwv, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")

ggplot(traits.RA, aes(y = chem.mass.mg_cwv, x = precip.cv)) +
  geom_smooth(method = "lm", se = T) + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Precip variability")


#### Time-Space: FD ####
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


## Functional richness ##
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


## Functional diversity
tmp <- as.data.frame(fd.all.ts$FDiv)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "FDiv"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp.sum <- tmp %>% 
  group_by(site, AI.site) %>%
  summarize(FDiv.mean = mean(FDiv, na.rm = T),
            FDiv.se = calcSE(FDiv))

ggplot(tmp.sum, aes(x = AI.site, y = FDiv.mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymax = FDiv.mean + FDiv.se, ymin = FDiv.mean - FDiv.se), width = 0.001)
 
## Functional dispersion
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

## Rao's Q
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

#### Fig 3 ####
fig3 <- ggarrange(fig3a, ggarrange(fig3b, fig3c, nrow = 1, ncol = 2, labels = c("(b)", "(c)")), labels = c("(a)","",""), nrow = 2, ncol = 1, heights = c(2.5,1))

ggsave("Manuscript/Aridity/Fig-3.png", fig3, height = 7, width = 6, dpi = 600, units = "in")

## CWM Traits ##
tmp <- as.data.frame(fd.all.ts$CWM)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))

tmp <- tmp %>%
  group_by(site) %>%
  summarize(across(chem.mass.mg:ldd.natural, ~ mean(.x, na.rm = T)))

tmp <- merge(tmp, meta[,c(1,2,19:22)], by.x = "site", by.y = "Dataset")

ggplot(tmp, aes(y = shape, x = log(AI.site))) +
  geom_smooth(method = "lm", se = T) + 
  #geom_point(size = 2) +
  labs(x = "Aridity Index") +
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(y = "CWM Shape")

summary(lm(shape~ log(AI.site), data = tmp))

ggplot(tmp, aes(y = size.mm, x = log(AI.site))) +
  geom_smooth(method = "lm", se = T) + 
  #geom_point(size = 2) +
  labs(x = "Aridity Index") +
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(y = "CWM Size")

summary(lm(shape~ log(AI.site), data = tmp))


### BY YEAR ####
## Functional richness ##
tmp <- as.data.frame(fd.all.ts$FRic)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "FRic"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp.sum <- tmp %>% 
  group_by(site, year, AI.site) %>%
  summarize(FRic.mean = mean(FRic, na.rm = T),
            FRic.se = calcSE(FRic))

ggplot(tmp.sum, aes(x = year, y = FRic.mean, col = site)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymax = FRic.mean + FRic.se, ymin = FRic.mean - FRic.se), width = 0.001)
 
m1 <- lm(FRic.mean ~ AI.site, data = tmp.sum)
summary(m1)

m1 <- lm(FRic ~ AI.site, data = tmp)
summary(m1) #YESSSS


## Functional diversity
tmp <- as.data.frame(fd.all.ts$FDiv)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "FDiv"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp.sum <- tmp %>% 
  group_by(site, year, AI.site) %>%
  summarize(FDiv.mean = mean(FDiv, na.rm = T),
            FDiv.se = calcSE(FDiv))

ggplot(tmp.sum, aes(x = year, y = FDiv.mean, col = site)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymax = FDiv.mean + FDiv.se, ymin = FDiv.mean - FDiv.se), width = 0.001)
 
## Functional dispersion
tmp <- as.data.frame(fd.all.ts$FDis)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "FDis"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp.sum <- tmp %>% 
  group_by(site, year, AI.site) %>%
  summarize(FDis.mean = mean(FDis, na.rm = T),
            FDis.se = calcSE(FDis))

ggplot(tmp.sum, aes(x = year, y = FDis.mean, col = site)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymax = FDis.mean + FDis.se, ymin = FDis.mean - FDis.se), width = 0.001)

## Rao's Q
tmp <- as.data.frame(fd.all.ts$RaoQ)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "RaoQ"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp.sum <- tmp %>% 
  group_by(site, year, AI.site) %>%
  summarize(RaoQ.mean = mean(RaoQ, na.rm = T),
            RaoQ.se = calcSE(RaoQ))

ggplot(tmp.sum, aes(x = year, y = RaoQ.mean, col = site)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymax = RaoQ.mean + RaoQ.se, ymin = RaoQ.mean - RaoQ.se), width = 0.001)

### BY space ####
## Functional richness ##
tmp <- as.data.frame(fd.all.ts$FRic)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "FRic"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp.sum <- tmp %>% 
  group_by(site, , AI.site) %>%
  summarize(FRic.mean = mean(FRic, na.rm = T),
            FRic.se = calcSE(FRic))

ggplot(tmp.sum, aes(x = year, y = FRic.mean, col = site)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymax = FRic.mean + FRic.se, ymin = FRic.mean - FRic.se), width = 0.001)
 
m1 <- lm(FRic.mean ~ AI.site, data = tmp.sum)
summary(m1)

m1 <- lm(FRic ~ AI.site, data = tmp)
summary(m1) #YESSSS


## Functional diversity
tmp <- as.data.frame(fd.all.ts$FDiv)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "FDiv"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp.sum <- tmp %>% 
  group_by(site, year, AI.site) %>%
  summarize(FDiv.mean = mean(FDiv, na.rm = T),
            FDiv.se = calcSE(FDiv))

ggplot(tmp.sum, aes(x = year, y = FDiv.mean, col = site)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymax = FDiv.mean + FDiv.se, ymin = FDiv.mean - FDiv.se), width = 0.001)
 
## Functional dispersion
tmp <- as.data.frame(fd.all.ts$FDis)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "FDis"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp.sum <- tmp %>% 
  group_by(site, year, AI.site) %>%
  summarize(FDis.mean = mean(FDis, na.rm = T),
            FDis.se = calcSE(FDis))

ggplot(tmp.sum, aes(x = year, y = FDis.mean, col = site)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymax = FDis.mean + FDis.se, ymin = FDis.mean - FDis.se), width = 0.001)

## Rao's Q
tmp <- as.data.frame(fd.all.ts$RaoQ)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "RaoQ"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp.sum <- tmp %>% 
  group_by(site, year, AI.site) %>%
  summarize(RaoQ.mean = mean(RaoQ, na.rm = T),
            RaoQ.se = calcSE(RaoQ))

ggplot(tmp.sum, aes(x = year, y = RaoQ.mean, col = site)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymax = RaoQ.mean + RaoQ.se, ymin = RaoQ.mean - RaoQ.se), width = 0.001)

#### OLD ####
#### Pairs Plots ####
source("/Users/Marina/Documents/UC-Davis/Quantitative-Resources/Scripts/Better_Pairs.R")
#hmm.trait <- filter(hmm.trait, site != "Carrizo")
colnames(hmm.trait)

# pairs(hmm.trait[hmm.trait$site == "Sonoran Desert", c(6,7, 21:23, 26:28)], 
#       lower.panel=panel.lmline, upper.panel=panel.cor)
# 
# pairs(hmm.trait[hmm.trait$site == "Sonoran Desert", c(6,7, 29,30, 32:34, 38)], 
#       lower.panel=panel.lmline, upper.panel=panel.cor)

pairs(hmm.trait[hmm.trait$system == "grassland", c(6,7, 21:23, 26:28)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

pairs(hmm.trait[hmm.trait$system == "grassland", c(6,7, 29,30, 32:34, 38)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)


pairs(hmm.trait[hmm.trait$system == "desert", c(6,7, 21:23, 26:28)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

pairs(hmm.trait[hmm.trait$system == "desert", c(6,7, 29, 30, 32:34, 38)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

## Grassland forb
### c: height (but... negative? higher things have less chance to go further? that makes NO sense)
### s: shape, size, c:n

## Desert forb
### c: shape, propC, 
### s: prop C, coat thick per width

## Grassland forb + grass
### c: shape, size, propC, E:S
### s: NOTHING

## Desert forb + grass
### c: shape, propC, E:S
### s: proportion C


#### AI Pairs ####
traits.forbs <- traits[traits$nat.inv == "native",]

# positive: mass, settling time, prop.N, wing loading, coat thickness, ldd
# species in less arid sites have heavier seeds that fall faster, higher nitrogen, higher coat thickness, disperse further

pairs(traits[, c(38, 40, 11:15)], 
      lower.panel = panel.lmline, 
      upper.panel = panel.cor)

pairs(traits[, c(38, 40, 16:20)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

pairs(traits[, c(38, 40, 21:25)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

pairs(traits[, c(38, 40, 26:27, 35, 36)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)



# ok but shouldnt community weighted traits matter more? just because they are the same doesnt mean shit

traits <- read.csv("Data/test.csv")

traits$AI.sum <- factor(traits$AI.sum, levels = c("arid", "arid-preferred", "semi-arid preferred", "semi-arid"))

traits$AI.bi <- ifelse(traits$AI.sum == "arid" | traits$AI.sum == "arid-preferred", "arid", "semi-arid")

traits.sum <- traits %>%
  filter(group == "forb") %>%
  group_by(nat.inv, AI.bi) %>%
  dplyr::summarize(across(morph.mass.mg:coat.thick.per.width,
                          list(mean = ~mean(.x, na.rm = TRUE),
                               se = ~ calcSE(.x))))


traits.sum <- traits %>%
  group_by(AI.sum) %>%
  dplyr::summarize(across(morph.mass.mg:coat.thick.per.width,
                          ~var(.x, na.rm = TRUE)))

ggplot(traits[traits$group == "forb",], aes(x = AI.sum, y = shape)) +
  geom_boxplot() +
  facet_wrap(~nat.inv)

ggplot(traits[traits$group == "forb",], aes(x = AI.sum, y = size.mm)) +
  geom_boxplot() +
  facet_wrap(~nat.inv)

ggplot(traits, aes(x = AI.sum, y = coat.thick.per.width)) +
  geom_boxplot() +
  facet_wrap(~nat.inv)

ggplot(traits, aes(x = AI.sum, y = prop.C)) +
  geom_boxplot() +
  facet_wrap(~nat.inv)

ggplot(traits, aes(x = AI.sum, y = prop.N)) +
  geom_boxplot() +
  facet_wrap(~nat.inv)

ggplot(traits, aes(x = AI.sum, y = set.time.mpsec)) +
  geom_boxplot() +
  facet_wrap(~nat.inv) 

ggplot(traits.sum, aes(x = AI.sum, y = shape)) +
  geom_point()

ggplot(traits.sum, aes(x = AI.sum, y = size.mm)) +
  geom_point() 

ggplot(traits.sum, aes(x = AI.sum, y = coat.thick.per.width)) +
  geom_point() 

ggplot(traits.sum, aes(x = AI.bi, y = coat.thick.per.width_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = coat.thick.per.width_mean - coat.thick.per.width_se, ymax = coat.thick.per.width_mean + coat.thick.per.width_se, width = 0.1)) +
  facet_wrap(~nat.inv)

ggplot(traits.sum, aes(x = AI.bi, y = both.thick_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = both.thick_mean - both.thick_se, ymax = both.thick_mean + both.thick_se, width = 0.1)) +
  facet_wrap(~nat.inv)

ggplot(traits.sum, aes(x = AI.bi, y = prop.C_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = prop.C_mean - prop.C_se, ymax = prop.C_mean + prop.C_se, width = 0.1)) +
  facet_wrap(~nat.inv)

ggplot(traits.sum, aes(x = AI.bi, y = prop.N_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = prop.N_mean - prop.N_se, ymax = prop.N_mean + prop.N_se, width = 0.1)) +
  facet_wrap(~nat.inv)

ggplot(traits.sum, aes(x = AI.bi, y = set.time.mpsec_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = set.time.mpsec_mean - set.time.mpsec_se, ymax = set.time.mpsec_mean + set.time.mpsec_se, width = 0.1)) +
  facet_wrap(~nat.inv)

ggplot(traits.sum, aes(x = AI.bi, y = ldd_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = ldd_mean - ldd_se, ymax = ldd_mean + ldd_se, width = 0.1)) +
  facet_wrap(~nat.inv)

ggplot(traits.sum, aes(x = AI.bi, y = morph.mass.mg_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = morph.mass.mg_mean - morph.mass.mg_se, ymax = morph.mass.mg_mean + morph.mass.mg_se, width = 0.1)) +
  facet_wrap(~nat.inv)

ggplot(traits.sum, aes(x = AI.sum, y = chem.mass.mg_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = chem.mass.mg_mean - chem.mass.mg_se, ymax = chem.mass.mg_mean + chem.mass.mg_se, width = 0.1)) +
  facet_wrap(~nat.inv)

# semi arid species are heavier, have thicker seed coats, 

# look at aridity index over time using 
# cv for precip or sd/var for temp
# mean and dispersion of traits and relate to mean and variability in aridity


