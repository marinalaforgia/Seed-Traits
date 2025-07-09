#### HMM System Comparison ####
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
library(rdacca.hp)
library(FD)
library(raster)
library(sf)
library(factoextra)
library(RColorBrewer)

#### Prep Long term Data ####
species <- read.csv("Data/20211001_Full-Species-List.csv")

# Carrizo #
CP.df <- readRDS("Scripts/HMMs/Carrizo-Plain/CP_HMM_none-control_ungrazed-GKR_2008.RDS")
#CP.df <- filter(CP.df, Code != "TRIGRA") # removing TRIGRA because... ?? 

# Jasper Ridge #
jr.df <- readRDS("Scripts/HMMs/Jasper-Ridge/jr_HMM_05m_control-avg.RDS") 

# Sonoran Desert #
sd.df <- readRDS("Scripts/HMMs/Sonoran-Desert/SD_HMM_2001.RDS")

# McLaughlin #
mcl.df <- readRDS("Scripts/HMMs/McLaughlin/20221212_mcl_HMM_quadrat.RDS")[,-14]
mcl.df <- filter(mcl.df, n.plots >= 80)

# Portal #
port.df <- readRDS("Scripts/HMMs/Portal/Port_HMM_control_winter_2009.RDS")

# skinner #
#skin.df <- readRDS("Scripts/HMMs/Skinner/skinner-hmm.RDS")

# Add sites #
sd.df$site <- "Sonoran"
port.df$site <- "Portal"
mcl.df$site <- "McLaughlin"
CP.df$site <- "Carrizo"
jr.df$site <- "Jasper Ridge"
#skin.df$site <- "Skinner"

meta <- read.csv("Data/Long-Term-Datasets/Datasets_metadata.csv")


# bind together #
full <- rbind(CP.df, jr.df, mcl.df, sd.df, port.df)

rm(CP.df, jr.df, mcl.df, sd.df, port.df)

full <- merge(full, meta[,c(1,2,19)], by.x = "site", by.y = "Dataset")

full <- filter(full, iter < 150)

#full <- filter(full, s >0)

#### Prep Seed Traits ####
traits <- read.csv("Data/20230801_Seed-Traits_clean_site.csv")

species.AI <- read.csv("Data/20230530_Seeds_All-Accessions.csv")

# AMSMEN has a suspiciously low prop.C compared to other Amsinckia species, replace with mean
traits[traits$Species == "Amsinckia menziesii",]$prop.C <- mean(traits$prop.C, na.rm = T)

traits <- filter(traits, new.code != "CUCPAL") # remove because ??
traits[is.na(traits$prop.C),]$prop.C <- mean(traits$prop.C, na.rm = T) # HERHIR
traits[is.na(traits$prop.N),]$prop.N <- mean(traits$prop.N, na.rm = T) # HERHIR
traits[is.na(traits$coat.perm.perc),]$coat.perm.perc <- mean(traits$coat.perm.perc, na.rm = T) # HERHIR

#traits <- filter(traits, new.code %in% full$Code)

# Normalize all traits
#hist(log(traits$wing.loading))
traits$wing.loading <- log(traits$wing.loading)

#hist(sqrt(traits$prop.C))
#traits$prop.C <- sqrt(traits$prop.C)

#hist(log(traits$prop.N))
#traits$prop.N <- log(traits$prop.N)

#hist(log(traits$coat.perm.perc))
traits$coat.perm.perc <- log(traits$coat.perm.perc)

#hist(log(traits$morph.mass.mg))
traits$morph.mass.mg <- log(traits$morph.mass.mg)

#hist(log(traits$chem.mass.mg))
traits$chem.mass.mg <- log(traits$chem.mass.mg)

#hist(log(traits$size.mm))
traits$size.mm <- log(traits$size.mm)

#hist(traits$set.time.mpsec)

#hist(sqrt(traits$shape))
#traits$shape <- sqrt(traits$shape)

#hist(log(traits$E.S))
#traits$E.S <- log(traits$E.S)

#hist(log(traits$coat.thick.per.width))
#traits$coat.thick.per.width <- log(traits$coat.thick.per.width)

#hist(traits$ldd.all)

#hist(traits$ldd.natural)

traits$FunGroup <- paste(traits$nat.inv, traits$group, sep = " ")

traits$FunGroup <- ifelse(traits$FunGroup == "native forb", "Native Forb", ifelse(traits$FunGroup == "native grass", "Native Grass", traits$FunGroup))

traits$FunGroup <- ifelse(traits$FunGroup == "invasive forb", "Exotic Forb", ifelse(traits$FunGroup == "invasive grass", "Exotic Grass", traits$FunGroup))

traits.ID <- traits

traits <- traits %>%
  group_by(Species, new.code, family, group, nat.inv, FunGroup, appendage.type, appendage, disp.cat.all, disp.cat.nat)  %>%
  dplyr::summarize(across(morph.mass.mg:carb.prop.kew, ~ mean(.x, na.rm = TRUE)))


# traits$AI.group <- ifelse(traits$AI.mean <= 0.2 & traits$AI.75 <= 0.2, "arid", 
#                           ifelse(traits$AI.mean > 0.2 & traits$AI.25 > 0.2, "semi-arid",
#                                  ifelse(traits$AI.mean > 0.2 & traits$AI.25 <= 0.2, "semi-arid preferred", "arid preferred")))
# 
# traits$AI.group <- factor(traits$AI.group, levels = c("arid", "arid preferred", "semi-arid preferred", "semi-arid"))

#traits <- traits[,c(1:26, 41, 27:40, 42)]

trait.col <- c("Species", "new.code", "family", "group", "nat.inv", "FunGroup", "morph.mass.mg", "chem.mass.mg",  "shape", "size.mm", "appendage.type", "appendage", "prop.C", "prop.N", "set.time.mpsec", "height.cm", "coat.perm.perc", "E.S", "both.thick", "mucilage", "fleshy.end", "starchy.end", "coat.thick.per.width", "ldd.all", "ldd.natural", "AI", "AI.mean", "AI.IQR", "AI.25", "AI.75", "AI.group")

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

k.cp$k.cat <- ifelse(k.cp$k.means == 1, "low S - high T",
                     ifelse(k.cp$k.means == 2, "high S - low T",
                            "medium S - medium T"))

### Portal ###
set.seed(123)
fviz_nbclust(hmm.trait[hmm.trait$site == "Portal",c(7,8)], kmeans, method = "wss") # 3

set.seed(123)
(k.portal <- kmeans(hmm.trait[hmm.trait$site == "Portal",c(7,8)], 3, nstart = 25))
# 1: high S - medium T
# 2: low S - high T
# 3: high S - low T

k.portal <- data.frame(Code = hmm.trait[hmm.trait$site == "Portal",]$Code, k.means = k.portal$cluster, site = "Portal")

k.portal$k.cat <- ifelse(k.portal$k.means == 1, "medium S - medium T",
                         ifelse(k.portal$k.means == 2, "low S - high T", 
                                "high S - low T"))

### Sonoran ###

set.seed(123)
fviz_nbclust(hmm.trait[hmm.trait$site == "Sonoran",c(7,8)], kmeans, method = "wss") # 3

set.seed(123)
(k.sd <- kmeans(hmm.trait[hmm.trait$site == "Sonoran",c(7,8)], 3, nstart = 25))

k.sd <- data.frame(Code = hmm.trait[hmm.trait$site == "Sonoran",]$Code, k.means = k.sd$cluster, site = "Sonoran")
# 1: medium S - high T
# 2: high S - low T
# 3: low S - high T
k.sd$k.cat <- ifelse(k.sd$k.means == 1, "medium S - medium T",
                     ifelse(k.sd$k.means == 2, "high S - low T",
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

### Skinner ###
# set.seed(123)
# fviz_nbclust(hmm.trait[hmm.trait$site == "Skinner",c(7,8)], kmeans, method = "wss") #3
# 
# set.seed(123)
# (k.skin <- kmeans(hmm.trait[hmm.trait$site == "Skinner",c(7,8)], 3, nstart = 25))
# # 1: medium S - medium T
# # 2: high S - low T
# # 3: low S - high T 
# k.skin <- data.frame(Code = hmm.trait[hmm.trait$site == "Skinner",]$Code, k.means = k.skin$cluster, site = "Skinner")
# 
# k.skin$k.cat <- ifelse(k.skin$k.means == 1, "medium S - medium T",
#                      ifelse(k.skin$k.means == 2, "low S - high T",
#                             "high S - low T"))

## bind together ##
k.means <- rbind(k.mcl, k.cp, k.jr, k.portal, k.sd)

hmm.trait <- merge(hmm.trait, k.means, by = c("Code", "site"))

#### Trade-off Graph ####
ggplot(hmm.trait, aes(x = s, y = c)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  theme_bw() +
  labs(x = "temporal dispersal", y = "spatial dispersal")

ggplot(hmm.trait, aes(x = s, y = c)) +
  geom_smooth(method = "lm", se = F, col = "black") + 
  geom_point(aes(col = k.cat, shape = k.cat), size = 2) +
  theme_bw() +
  theme(
    legend.position = c(0.83, 0.3),
    legend.title = element_blank()
  ) +
  facet_wrap(~site, scales = "free") +
  scale_color_viridis_d() +
  labs(x = "temporal dispersal", y = "spatial dispersal")


k.cat <- hmm.trait %>%
  dplyr::group_by(site) %>%
  dplyr::mutate(total = length(Species)) %>%
  dplyr::group_by(site, k.cat, total) %>%
  dplyr::summarize(by.cat = length(Species)) %>%
  dplyr::mutate(cat.pro = by.cat/total)

k.cat$k.cat <- factor(k.cat$k.cat, levels = c("low S - high T", "medium S - medium T", "high S - low T"))

ggplot(k.cat, aes(x = site, y = cat.pro, group = k.cat, fill = k.cat)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(y = "Relative Proportion of Community")

#### PCA ####
#### .___PCA ALL ####
trait.col <- c("Species", "Code", "family", "group", "nat.inv", "FunGroup", "morph.mass.mg", "chem.mass.mg",  "shape", "size.mm", "appendage.type", "appendage", "prop.C", "prop.N", "set.time.mpsec", "height.cm", "coat.perm.perc", "E.S", "both.thick", "mucilage", "fleshy.end", "starchy.end", "coat.thick.per.width", "wing.loading", "ldd.all", "ldd.natural", "AI", "AI.mean", "AI.IQR", "AI.25", "AI.75", "AI.group", "k.cat")

trait.pca.col <- c("chem.mass.mg",  "shape", 
"prop.C", "prop.N", "set.time.mpsec", "E.S", "size.mm", "ldd.natural", "coat.thick.per.width", "height.cm", "wing.loading")

trait.pca <- hmm.trait[complete.cases(hmm.trait[,"Species"]), trait.col]
pca <- prcomp(trait.pca[, trait.pca.col], scale = T)
summary(pca)

hmm.trait <- cbind(hmm.trait[complete.cases(hmm.trait[,"Species"]),], pca$x[,1:4])

#trait.pca[is.na(trait.pca$site),]$site <- "none"

# trait.pca$AI.group2 <- ifelse(trait.pca$AI.group == "arid" | trait.pca$AI.group == "arid preferred", "arid", "semi-arid")
# 
#trait.pca$site2 <- ifelse(trait.pca$site == "Jasper Ridge" | trait.pca$site == "McLaughlin", "semi-arid", "arid")
# 
# trait.pca$site2 <- ifelse(trait.pca$site == "Portal" | trait.pca$site == "Sonoran", "arid", "semi-arid")

autoplot(pca, x = 1, y = 2, data = hmm.trait, frame = F, loadings = T, loadings.label = T, label = F, col = "k.cat", size = 3) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
    ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"), na.translate = F) +
  stat_ellipse(aes(group = k.cat, col = k.cat)) 
  #facet_wrap(~site2) #+
  #geom_text(data = trait.pca[!is.na(trait.pca$site),], aes(label = Code))



m1 <- lm(PC2 ~ site, trait.pca[trait.pca$site != "none",])
anova(m1)

autoplot(pca, x = 1, y = 2, data = hmm.trait, frame = F, loadings = T, loadings.label = T, label = F, col = "site", size = 2) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"), na.translate = F) +
  stat_ellipse(aes(group = site, col = site))

autoplot(pca, x = 3, y = 4, data = hmm.trait, frame = F, loadings = T, loadings.label = T, label = F, col = "site", size = 2) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"), na.translate = F) +
  stat_ellipse(aes(group = site, col = site))

pca_hull <- 
  trait.pca %>% 
  filter(site != "none") %>%
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

ggplot(trait.pca[trait.pca$site != "none",], aes(x = PC1, y = PC2)) +
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


### SAC STATE GRAPHS ##
# talk <- autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = F, label = F, col = "darkgrey", size = 0.5) +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(colour = "black", fill = NA, size = 1),
#     legend.title = element_text(size = 15),
#     legend.text = element_text(size = 14),
#     axis.title = element_text(size = 15),
#     axis.text = element_text(size = 14)
#   ) #+
#   # scale_color_manual(values = c(adhesive = "#1B9E77", ant = "red4", unassisted = "darkgoldenrod3", ingestion = "#7570B3", wind = "#F17236"))
# 
# talk.ant <- autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = F, label = F, col = "ldd2", size = 2.5, shape = "ldd2", alpha = "ldd2") +
#   #geom_point(aes(size = 2)) +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(colour = "black", fill = NA, size = 1),
#     legend.title = element_text(size = 15),
#     legend.text = element_text(size = 14),
#     axis.title = element_text(size = 15),
#     axis.text = element_text(size = 14)
#   ) +
#   scale_color_manual(values = c(adhesive = "#1B9E77", ant = "red4", unassisted = "darkgoldenrod3", ingestion = "#7570B3", wind = "#F17236")) +
#   scale_alpha_manual(values = c(adhesive = 0, ant = 1, unassisted = 0, ingestion = 0, wind = 0))
# 
# talk.adhesive <- autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = F, label = F, col = "ldd2", size = 2.5, shape = "ldd2", alpha = "ldd2") +
#   #geom_point(aes(size = 2)) +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(colour = "black", fill = NA, size = 1),
#     legend.title = element_text(size = 15),
#     legend.text = element_text(size = 14),
#     axis.title = element_text(size = 15),
#     axis.text = element_text(size = 14)
#   ) +
#   scale_color_manual(values = c(adhesive = "#1B9E77", ant = "red4", unassisted = "darkgoldenrod3", ingestion = "#7570B3", wind = "#F17236")) +
#   scale_alpha_manual(values = c(adhesive = 1, ant = 0, unassisted = 0, ingestion = 0, wind = 0))
# 
# talk.wind <- autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = F, label = F, col = "ldd2", size = 2.5, shape = "ldd2", alpha = "ldd2") +
#   #geom_point(aes(size = 2)) +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(colour = "black", fill = NA, size = 1),
#     legend.title = element_text(size = 15),
#     legend.text = element_text(size = 14),
#     axis.title = element_text(size = 15),
#     axis.text = element_text(size = 14)
#   ) +
#   scale_color_manual(values = c(adhesive = "#1B9E77", ant = "red4", unassisted = "darkgoldenrod3", ingestion = "#7570B3", wind = "#F17236")) +
#   scale_alpha_manual(values = c(adhesive = 0, ant = 0, unassisted = 0, ingestion = 0, wind = 1))
# 
# talk.un <- autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = F, label = F, col = "ldd2", size = 2.5, shape = "ldd2", alpha = "ldd2") +
#   #geom_point(aes(size = 2)) +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(colour = "black", fill = NA, size = 1),
#     legend.title = element_text(size = 15),
#     legend.text = element_text(size = 14),
#     axis.title = element_text(size = 15),
#     axis.text = element_text(size = 14)
#   ) +
#   scale_color_manual(values = c(adhesive = "#1B9E77", ant = "red4", unassisted = "darkgoldenrod3", ingestion = "#7570B3", wind = "#F17236")) +
#   scale_alpha_manual(values = c(adhesive = 0, ant = 0, unassisted = 1, ingestion = 0, wind = 0))
# 
# talk.ing <- autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = F, label = F, col = "ldd2", size = 2.5, shape = "ldd2", alpha = "ldd2") +
#   #geom_point(aes(size = 2)) +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(colour = "black", fill = NA, size = 1),
#     legend.title = element_text(size = 15),
#     legend.text = element_text(size = 14),
#     axis.title = element_text(size = 15),
#     axis.text = element_text(size = 14)
#   ) +
#   scale_color_manual(values = c(adhesive = "#1B9E77", ant = "red4", unassisted = "darkgoldenrod3", ingestion = "#7570B3", wind = "#F17236")) +
#   scale_alpha_manual(values = c(adhesive = 0, ant = 0, unassisted = 0, ingestion = 1, wind = 0))
# 
# talk.blank <- autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = F, label = F, col = "white", size = 2.5, shape = "ldd2") +
#   #geom_point(aes(size = 2)) +
#   theme_classic() +
#   theme(
#     panel.border = element_rect(colour = "black", fill = NA, size = 1),
#     legend.title = element_text(size = 15),
#     legend.text = element_text(size = 14),
#     axis.title = element_text(size = 15),
#     axis.text = element_text(size = 14)
#   )
# 
# ggsave("Extended-Figures/talk.jpeg", talk, width = 4.5, height = 4)
# 
# ggsave("Extended-Figures/talk-blank.jpeg", talk.blank, width = 6, height = 4)
# 
# ggsave("Extended-Figures/talk-adhesive.jpeg", talk.adhesive, width = 6, height = 4)
# 
# ggsave("Extended-Figures/talk-ant.jpeg", talk.ant, width = 6, height = 4)
# 
# ggsave("Extended-Figures/talk-ing.jpeg", talk.ing, width = 6, height = 4)
# 
# ggsave("Extended-Figures/talk-un.jpeg", talk.un, width = 6, height = 4)
# 
# ggsave("Extended-Figures/talk-wind.jpeg", talk.wind, width = 6, height = 4)
# 
# library(Rmisc)
# trait.pca.sum <- summarySE(trait.pca, groupvars = "ldd2", measurevar = "PC2")
# 
# 
# ggplot(trait.pca.sum[trait.pca.sum$N >3,], aes(x = ldd2, y = PC2)) +
#   geom_point(size = 2) +
#   geom_errorbar(aes(ymin = PC2 - se, ymax = PC2 + se), width = 0.2) +
#   theme_classic() +
#   theme(
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 15)
#   )
# 
# trait.pca <- merge(trait.pca, disp.sum[,c(1,3)], by = "Species", all.x = T)
# 
# trait.pca.sum <- summarySE(trait.pca, groupvars = "ldd", measurevar = "PC2")
# 
# disp.fig <- ggplot(trait.pca.sum[complete.cases(trait.pca.sum),], aes(x = ldd, y = PC2)) +
#   geom_point(size = 2) +
#   geom_errorbar(aes(ymin = PC2 - se, ymax = PC2 + se), width = 0.1) +
#   theme_classic() +
#   theme(
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 15)
#   )
# 
# ggsave("Extended-Figures/disp.png", disp.fig , units = "in", height = 3, width = 5.5)
#### END

#### .___PCA Forbs ####
trait.pca.col <- c("chem.mass.mg",  "shape", "prop.C", "prop.N", "set.time.mpsec", "coat.perm.perc", "E.S", "size.mm", "ldd.natural", "coat.thick.per.width", "height.cm")

hmm.trait.forb <- merge(hmm.trait.forb, hmm.trait[,c(1,58)], by = "Code", all.x = T, all.y = F)
trait.pca <- hmm.trait.forb[complete.cases(hmm.trait.forb[,"Species"]), trait.col]
pca <- prcomp(trait.pca[, trait.pca.col], scale = T)
summary(pca)

hmm.trait.forb <- cbind(hmm.trait.forb[complete.cases(hmm.trait.forb[,"Species"]),], pca$x[,1:4])
#trait.pca <- merge(trait.pca, hmm.trait[,c(1,2,53)], by = c("Code", "site"), all.x = T, all.y = F)
#trait.pca[is.na(trait.pca$site),]$site <- "none"

# trait.pca$AI.group2 <- ifelse(trait.pca$AI.group == "arid" | trait.pca$AI.group == "arid preferred", "arid", "semi-arid")
# 
#trait.pca$site2 <- ifelse(trait.pca$site == "Jasper Ridge" | trait.pca$site == "McLaughlin", "semi-arid", "arid")
# 
# trait.pca$site2 <- ifelse(trait.pca$site == "Portal" | trait.pca$site == "Sonoran", "arid", "semi-arid")

autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "k.cat", size = 3) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
    ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"), na.translate = F) +
  stat_ellipse(aes(group = k.cat, col = k.cat)) 
  #facet_wrap(~site2) #+
  #geom_text(data = trait.pca[!is.na(trait.pca$site),], aes(label = Code))

autoplot(pca, x = 3, y = 4, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "k.cat", size = 2) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"), na.translate = F) +
  #scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"), na.translate = F) +
  stat_ellipse(aes(group = k.cat, col = k.cat))

ggplot(trait.pca.forb, aes(x = ldd.natural, y = shape)) +
  geom_point() +
  geom_smooth(method = "lm")

####.___ JUST species in DF ####
trait.pca.col <- c("chem.mass.mg",  "shape", 
"prop.C", "prop.N", "set.time.mpsec", "coat.perm.perc", "E.S", "size.mm", "ldd.natural", "coat.thick.per.width")

traits.site <- merge(traits, full[full$site != "Skinner", c(1,3,6,7,8,15,16)], by.x = "new.code", by.y = "Code")

#traits.site <- traits.site[traits.site$site != "Skinner", c("site", trait.col)]

#traits.site <- traits.site[traits.site$FunGroup == "Native Forb",]
pca <- prcomp(traits.site[, trait.pca.col], scale = T)
summary(pca)

traits.site <- cbind(traits.site, pca$x[,1:4])

#brewer.pal(n = 5, name = "Dark2")

autoplot(pca, x = 1, y = 2, data = traits.site, frame = F, loadings = T, loadings.label = T, label = F, size = 2, col = "site", shape =) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
    )  +
  stat_ellipse(aes(group = site, col = site))#+
   #geom_text(aes(label = new.code))

autoplot(pca, x = 3, y = 4, data = traits.site, frame = F, loadings = T, loadings.label = T, label = F, size = 2, col = "site", shape = "site") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
    )

m1 <- lm(PC1 ~ site, data = traits.site)
#plot(m1)
anova(m1)

m1 <- lm(PC2 ~ site, data = traits.site)
#plot(m1)
anova(m1)

m1 <- lm(PC3 ~ site, data = traits.site)
#plot(m1)
anova(m1)

m1 <- lm(PC4 ~ site, data = traits.site)
#plot(m1)
anova(m1)

# PC1 is grasses v forbs
# PC2 is aridity gradient

pca_hull <- 
  traits.site %>% 
  group_by(group) %>% 
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

ggplot(traits.site, aes(x = PC1, y = PC2)) +
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

pca_hull <- 
  traits.site %>% 
  group_by(group) %>% 
  slice(chull(PC1, PC2))

ggplot(traits.site, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = group)) +
  theme_light() + 
  geom_polygon(data = pca_hull,
               aes(fill = group,
                   colour = group),
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

autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, size = 2, col = "FunGroup", shape = "FunGroup") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
    )

autoplot(pca, x = 3, y = 4, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, size = 2, col = "site", shape = "site") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
    )

#autoplot(kmeans(trait.pca[,trait.pca.col], 3), data = trait.pca[,trait.pca.col])

####.___ JUST Forbs in DF ####
trait.pca.col <- c("chem.mass.mg",  "shape", 
"prop.C", "prop.N", "set.time.mpsec", "coat.perm.perc", "E.S", "size.mm", "ldd.natural", "coat.thick.per.width")

traits.site <- merge(traits, full[full$site != "Skinner", c(1,3,6,7,8,15,16)], by.x = "new.code", by.y = "Code")


traits.site.f <- traits.site[traits.site$group == "forb",]
pca <- prcomp(traits.site.f[, trait.pca.col], scale = T)
summary(pca)

traits.site.f <- cbind(traits.site.f, pca$x[,1:4])

#brewer.pal(n = 5, name = "Dark2")

autoplot(pca, x = 1, y = 2, data = traits.site.f, frame = F, loadings = T, loadings.label = T, label = F, size = 2, col = "site", shape =) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
    )  +
  stat_ellipse(aes(group = site, col = site))#+
   #geom_text(aes(label = new.code))

autoplot(pca, x = 3, y = 4, data = traits.site, frame = F, loadings = T, loadings.label = T, label = F, size = 2, col = "site", shape = "site") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
    )


###.___ PCA Arid ####
trait.pca.col <- c("chem.mass.mg",  "shape", 
"prop.C", "prop.N", "set.time.mpsec", "coat.perm.perc", "E.S", "size.mm", "coat.thick.per.width")

trait.pca <- traits[, trait.col]
trait.pca <- merge(full[,c(1,3)], trait.pca, by.y = "new.code" , by.x = "Code", all.x = F, all.y = T)

# trait.pca <- filter(trait.pca, site != "Sonoran", site != "Portal", group == "forb")
# 
# trait.pca <- filter(trait.pca, site == "Sonoran" | site == "Portal", group == "forb")

#trait.pca <- filter(trait.pca, AI.group == "arid", group == "forb")

#trait.pca <- filter(trait.pca, AI.group == "semi-arid", group == "forb")

pca <- prcomp(trait.pca[, trait.pca.col], scale = T)
summary(pca)

trait.pca <- cbind(trait.pca, pca$x[,1:4])
#trait.pca[is.na(trait.pca$site),]$site <- "none"

autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "site", size = 3) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
    ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"), na.translate = F) +
  stat_ellipse(aes(group = site, col = site)) 
 # geom_text(aes(label = Code))

#### .___PCA kmeans ####
trait.pca.col <- c("chem.mass.mg",  "shape", 
"prop.C", "prop.N", "set.time.mpsec", "coat.perm.perc", "E.S", "size.mm", "ldd.natural", "coat.thick.per.width", "height.cm")

#trait.pca <- traits[, trait.col]
#trait.pca <- merge(trait.pca, hmm.trait, by.x = "new.code", by.y = "Code")
pca <- prcomp(hmm.trait[, trait.pca.col], scale = T)
summary(pca)

hmm.trait <- cbind(hmm.trait[complete.cases(hmm.trait),-c(46:49)], pca$x[,1:4])
#trait.pca[is.na(trait.pca$site),]$site <- "none"

autoplot(pca, x = 1, y = 2, data = hmm.trait, frame = F, loadings = T, loadings.label = T, label = F, col = "k.cat", size = 3) +
  #geom_jitter() + 
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
    ) #+
  #scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"), na.translate = F) #+
  #stat_ellipse(aes(group = site, col = site))

#### PC graphs #### 

####.___PC1 ####

##FULL##
## mostly size (grasses) which explains why this pattern most holds true for mclaughlin and a little for JR
#hmm.trait <- traits.site
ggplot(hmm.trait, aes(x = PC1, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

# similarly here, mostly mclauhglin and JR, some sonoran
ggplot(hmm.trait, aes(x = PC1, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

##FORBS ONLY##
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

####.___PC2 ####
# no clear trait relationships, lots of traits (mass, settling time, coat perm, coat thickness)
##FULL##
ggplot(hmm.trait, aes(x = PC2, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

ggplot(hmm.trait, aes(x = PC2, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

##FORBS ONLY##
# PC2 IS MOSTLY SHAPE and dispersal potential, in most environments, rounder seeds have higher temporal dispersal, but not in the sonoran desert
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
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

# ggplot(traits.site.f, aes(x = PC2, y = c)) +
#   geom_smooth(method = "lm") + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "spatial dispersal")
# 
# ggplot(traits.site.f, aes(x = PC2, y = s)) +
#   geom_smooth(method = "lm") + 
#   #geom_text(aes(label = Code)) + 
#   geom_point() +
#   theme_bw() +
#   facet_wrap(~site, scales = "free") +
#   labs(y = "temporal dispersal")

hmm.trait.forb.sum <- hmm.trait.forb %>%
  group_by(site, AI.site) %>%
  summarize(PC2.mean = mean(PC2, na.rm = T),
            PC2.se = calcSE(PC2))
  
ggplot(hmm.trait.forb.sum, aes(x = log(AI.site), y = PC2.mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = PC2.mean - PC2.se, ymax = PC2.mean + PC2.se), width = 0.1)

####.___PC3####

##FULL##
ggplot(hmm.trait, aes(x = PC3, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

ggplot(hmm.trait, aes(x = PC3, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

##FORBS ONLY##
ggplot(hmm.trait.forb, aes(x = PC3, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

ggplot(hmm.trait.forb, aes(x = PC3, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

####.___PC4####

##FULL##
ggplot(hmm.trait, aes(x = PC4, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

ggplot(hmm.trait, aes(x = PC4, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

##FORBS ONLY##
ggplot(hmm.trait.forb, aes(x = PC4, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal")

ggplot(hmm.trait.forb, aes(x = PC4, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")



#### Hier Partition ####
survival <- c("shape", "prop.C", "prop.N", "coat.perm.perc", "coat.thick.per.width", "size.mm", "chem.mass.mg", "E.S")

colonization <- c("shape", "chem.mass.mg", "set.time.mpsec", "size.mm", "ldd.natural")

##### ALL ####
hmm.trait <- hmm.trait[complete.cases(hmm.trait[,-(33:35)]),]

# combining them all together leads to like... nothing
HPA.s <- rdacca.hp(
  dv = hmm.trait$s,
  iv = hmm.trait[, colnames(hmm.trait) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation #5% explained variation

plot.rdaccahp(HPA.s) # shape, size, proportion N

HPA.c <- rdacca.hp(
  dv = hmm.trait$c,
  iv = hmm.trait[, colnames(hmm.trait) %in% colonization],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.c$Hier.part
#HPA.s$Var.part
HPA.c$Total_explained_variation #5% explained variation


plot.rdaccahp(HPA.c) # ldd.natural, settling time 

###### by site ####
HPA.s <- rdacca.hp(
  dv = hmm.trait[hmm.trait$site == "Carrizo",]$s,
  iv = hmm.trait[hmm.trait$site == "Carrizo", colnames(hmm.trait) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation

plot.rdaccahp(HPA.s)  #prop C is the only positive one

HPA.s <- rdacca.hp(
  dv = hmm.trait[hmm.trait$site == "Sonoran",]$s,
  iv = hmm.trait[hmm.trait$site == "Sonoran", colnames(hmm.trait) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation

plot.rdaccahp(HPA.s) # shape and coat thickness

HPA.s <- rdacca.hp(
  dv = hmm.trait[hmm.trait$site == "Portal",]$s,
  iv = hmm.trait[hmm.trait$site == "Portal", colnames(hmm.trait) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation

plot.rdaccahp(HPA.s) # shape and coat thcikness

HPA.s <- rdacca.hp(
  dv = hmm.trait[hmm.trait$site == "Jasper Ridge",]$s,
  iv = hmm.trait[hmm.trait$site == "Jasper Ridge", colnames(hmm.trait) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation

plot.rdaccahp(HPA.s) # size and prop C

HPA.s <- rdacca.hp(
  dv = hmm.trait[hmm.trait$site == "McLaughlin",]$s,
  iv = hmm.trait[hmm.trait$site == "McLaughlin", colnames(hmm.trait) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation

plot.rdaccahp(HPA.s) # shape and size

# Carrizo has neg explained amount, the rest have in the 50s (mc 15)
HPA.c <- rdacca.hp(
  dv = hmm.trait[hmm.trait$site == "Carrizo",]$c,
  iv = hmm.trait[hmm.trait$site == "Carrizo", colnames(hmm.trait) %in% colonization],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.c$Hier.part
#HPA.c$Var.part
HPA.c$Total_explained_variation

plot.rdaccahp(HPA.c) 

HPA.c <- rdacca.hp(
  dv = hmm.trait[hmm.trait$site == "Sonoran",]$c,
  iv = hmm.trait[hmm.trait$site == "Sonoran", colnames(hmm.trait) %in% colonization],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.c$Hier.part
#HPA.c$Var.part
HPA.c$Total_explained_variation

plot.rdaccahp(HPA.c) 

HPA.c <- rdacca.hp(
  dv = hmm.trait[hmm.trait$site == "Portal",]$c,
  iv = hmm.trait[hmm.trait$site == "Portal", colnames(hmm.trait) %in% colonization],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.c$Hier.part
#HPA.c$Var.part
HPA.c$Total_explained_variation

plot.rdaccahp(HPA.c) 

HPA.c <- rdacca.hp(
  dv = hmm.trait[hmm.trait$site == "Jasper Ridge",]$c,
  iv = hmm.trait[hmm.trait$site == "Jasper Ridge", colnames(hmm.trait) %in% colonization],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.c$Hier.part
#HPA.c$Var.part
HPA.c$Total_explained_variation

plot.rdaccahp(HPA.c) 

HPA.c <- rdacca.hp(
  dv = hmm.trait[hmm.trait$site == "McLaughlin",]$c,
  iv = hmm.trait[hmm.trait$site == "McLaughlin", colnames(hmm.trait) %in% colonization],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.c$Hier.part
#HPA.c$Var.part
HPA.c$Total_explained_variation

plot.rdaccahp(HPA.c) 

##### FORBS ####

hmm.trait.forb <- hmm.trait.forb[complete.cases(hmm.trait.forb[,-(33:35)]),]

HPA.s <- rdacca.hp(
  dv = hmm.trait.forb$s,
  iv = hmm.trait.forb[, colnames(hmm.trait.forb) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation

plot.rdaccahp(HPA.s) 

HPA.c <- rdacca.hp(
  dv = hmm.trait.forb$c,
  iv = hmm.trait.forb[, colnames(hmm.trait.forb) %in% colonization],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.c$Hier.part
#HPA.s$Var.part
HPA.c$Total_explained_variation

plot.rdaccahp(HPA.c) 

# it's so variable that nothing is really doing a good job of explaining dispersal through space or time

###### by site ####
hmm.trait.forb <- hmm.trait.forb[complete.cases(hmm.trait.forb[,-(33:35)]),]

HPA.s <- rdacca.hp(
  dv = hmm.trait.forb[hmm.trait.forb$site == "Carrizo",]$s,
  iv = hmm.trait.forb[hmm.trait.forb$site == "Carrizo", colnames(hmm.trait.forb) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation # 5% looool 

plot.rdaccahp(HPA.s) 

HPA.s <- rdacca.hp(
  dv = hmm.trait.forb[hmm.trait.forb$site == "Sonoran",]$s,
  iv = hmm.trait.forb[hmm.trait.forb$site == "Sonoran", colnames(hmm.trait.forb) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation #8%

plot.rdaccahp(HPA.s) 

HPA.s <- rdacca.hp(
  dv = hmm.trait.forb[hmm.trait.forb$site == "Portal",]$s,
  iv = hmm.trait.forb[hmm.trait.forb$site == "Portal", colnames(hmm.trait.forb) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation

plot.rdaccahp(HPA.s) 

HPA.s <- rdacca.hp(
  dv = hmm.trait.forb[hmm.trait.forb$site == "Jasper Ridge",]$s,
  iv = hmm.trait.forb[hmm.trait.forb$site == "Jasper Ridge", colnames(hmm.trait.forb) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation

plot.rdaccahp(HPA.s) 

HPA.s <- rdacca.hp(
  dv = hmm.trait.forb[hmm.trait.forb$site == "McLaughlin",]$s,
  iv = hmm.trait.forb[hmm.trait.forb$site == "McLaughlin", colnames(hmm.trait.forb) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation

plot.rdaccahp(HPA.s) 


#### Trait graphs #### 
## s: shape, size, prop.N
## c: settling time, ldd.natural

##### Spatial Dispersal ####

ggplot(hmm.trait, aes(x = ldd.natural, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "dispersal potential")

ggplot(hmm.trait, aes(x = ldd.all, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "dispersal potential")

ggplot(hmm.trait, aes(x = height.cm, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "height") # is height more important in SD??

ggplot(hmm.trait, aes(x = set.time.mpsec, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "settling time (mpsec)")

ggplot(hmm.trait, aes(x = chem.mass.mg, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "mass")

ggplot(hmm.trait, aes(x = morph.mass.mg, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "mass")

ggplot(hmm.trait, aes(x = size.mm, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "spatial dispersal", x = "mass")


##### Temporal Dispersal ####
ggplot(hmm.trait, aes(x = shape, y = s)) +
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


### OTHERS
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

ggplot(hmm.trait.forb, aes(x = shape, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  labs(y = "temporal dispersal")

ggplot(hmm.trait.forb, aes(x = size.mm, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
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

#### Trait CV ####

disp <- hmm.trait %>%
  group_by(site) %>%
  dplyr::summarize(across(chem.mass.mg:ldd.natural, ~sd(.x, na.rm = TRUE)/mean(.x, na.rm = T)))

disp <- merge(disp, meta[,c(1,2,19)], by.x = "site", by.y = "Dataset")

ggplot(disp, aes(x = precip.cv, y = coat.thick.per.width)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "precip cv", y = "coat thickness cv") 


ggplot(disp, aes(x = precip.cv, y = shape)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "precip cv", y = "shape cv") 

ggplot(disp, aes(x = precip.cv, y = size.mm)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "precip cv", y = "size cv") 
  
ggplot(disp, aes(x = precip.cv, y = chem.mass.mg)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "precip cv", y = "mass cv") 

disp$AI.site <- disp$AI.site*0.0001

ggplot(disp, aes(x = AI.site, y = chem.mass.mg)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "AI site", y = "mass cv") 

ggplot(disp, aes(x = AI.site, y = size.mm)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "AI site", y = "size cv") 

ggplot(disp, aes(x = AI.site, y = shape)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "AI site", y = "shape cv") 

ggplot(disp, aes(x = AI.site, y = coat.thick.per.width)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "AI site", y = "coat thick cv") 


ggplot(disp, aes(x = AI.site, y = prop.N)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "AI site", y = "N cv") 


ggplot(disp, aes(x = AI.site, y = prop.C)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "AI site", y = "N cv") 

ggplot(disp, aes(x = AI.site, y = set.time.mpsec)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "AI site", y = "settling time cv") 


ggplot(disp, aes(x = AI.site, y = ldd.natural)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "AI site", y = "dispersal potential cv") 

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

#### Trait dispersion across 300 species ####
traits.disp <- traits %>%
  filter(nat.inv == "native") %>%
  group_by(AI.group) %>%
  dplyr::summarize(across(morph.mass.mg:coat.thick.per.width, ~var(.x, na.rm = TRUE)))

count(traits$AI.group)


ggplot(traits.disp, aes(x = AI.group, y = chem.mass.mg)) +
  geom_point()

ggplot(traits.disp, aes(x = AI.group, y = shape)) +
  geom_point()

ggplot(traits.disp, aes(x = AI.group, y = ldd.natural)) +
  geom_point()

ggplot(traits.disp, aes(x = AI.group, y = coat.thick.per.width)) +
  geom_point()

ggplot(traits.disp, aes(x = AI.group, y = size.mm)) +
  geom_point()

ggplot(traits.disp, aes(x = AI.group, y = prop.N)) +
  geom_point()

ggplot(traits.disp, aes(x = AI.group, y = prop.C)) +
  geom_point()

ggplot(traits.disp, aes(x = AI.group, y = set.time.mpsec)) +
  geom_point()


#### CWM traits ####
mcl <- readRDS("Data/Long-Term-Datasets/McLaughlin/McL_RA.RDS")
mcl$site <- "McLaughlin"
names(mcl)

jr <- readRDS("Data/Long-Term-Datasets/Jasper-Ridge/JR_RA.rds")
jr$site <- "Jasper Ridge"

#skin <- readRDS("Data/Long-Term-Datasets/Skinner/skin-ra.RDS")
#skin$site <- "Skinner"

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
	

traits.RA <- merge(traits, df.RA, by.x = "new.code", by.y = "Code", all = F)

traits.RA <- traits.RA[,c(1:27,35,36,43,44)]

traits.RA <- merge(traits.RA, hmm.trait[,c(1,2,59:62)], by.y = c("Code", "site"), by.x = c("new.code", "site"))

# CWM and CWV (see Holshof et al 2013)

traits.RA <- traits.RA %>%
  group_by(site) %>%
  summarize(across(morph.mass.mg:PC4, 
                   list(cwm = ~ sum(.x * RA, na.rm = T),
                        cwv = ~ sum(RA * I(.x - sum(.x * RA, na.rm = T))^2, na.rm = T)
                          )))


traits.RA <- merge(traits.RA, meta[,c(1,2,19:22)], by.x = "site", by.y = "Dataset")

##### CWM: Aridity index ####

ggplot(traits.RA, aes(y = PC2_cwm, x = AI.site)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index") 

cor.test(traits.RA$PC2_cwm, traits.RA$AI.site)


ggplot(traits.RA, aes(y = PC2_cwv, x = AI.site)) +
  geom_smooth(method = "lm") + 
  geom_text(aes(label = site)) +
  theme_bw() +
  labs(x = "Aridity Index") 

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

b <- ggplot(traits.RA, aes(y = shape_cwm, x = AI.site)) +
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

#### NMDS ####

# other way of doing NMDS (perhaps more correct)
traits.RA.wide <- pivot_wider(traits.RA, names_from = Species, values_from = RA)
traits.RA.wide[is.na(traits.RA.wide)] <- 0
traits.RA.wide <- merge(traits.RA.wide, Plot.cwm[,c(1:7)], by = "Plot")

traits.RA.nmds <- metaMDS(traits.RA.wide[,3:94], k = 3, trymax = 100) 
stressplot(traits.RA.nmds)

# JPEG device
pdf("nmds_plot.pdf", width = 6, height = 5)

ordiplot(SB.sum.nmds)

#ordiellipse(SB.sum.nmds, groups = SB.sum.wide$Serpentine, draw = "polygon", col = c("red", "green", "blue"), label = T)

ordihull(SB.sum.nmds, groups = SB.sum.wide$Serpentine, draw = "polygon", col = c("red", "green", "blue"), label = T)

env <- SB.sum.wide[,c(95:100)]
colnames(env) <- c("Shape", "ST", "Mass", "Wing Loading", "Height", "LDD")
en <- envfit(SB.sum.nmds, env, permutations = 9999)
plot(en)

# Close device
dev.off()


#### Mean: FD ####
mcl <- readRDS("Data/Long-Term-Datasets/McLaughlin/McL_RA.RDS")
mcl$site <- "McLaughlin"
names(mcl)

jr <- readRDS("Data/Long-Term-Datasets/Jasper-Ridge/JR_RA.rds")
jr$site <- "Jasper Ridge"

#skin <- readRDS("Data/Long-Term-Datasets/Skinner/skin-ra.RDS")
#skin$site <- "Skinner"

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

df.RA.wide <- pivot_wider(df.RA, names_from = "Code", values_from = RA)
df.RA.wide <- as.data.frame(df.RA.wide)
row.names(df.RA.wide) <- df.RA.wide$site
df.RA.wide <- df.RA.wide[,-1]
df.RA.wide <- df.RA.wide[,order(colnames(df.RA.wide))]

traits.FD <- traits[,c(2,12:20,24,25,27)]
traits.FD <- as.data.frame(traits.FD)
row.names(traits.FD) <- traits.FD$new.code
traits.FD <- traits.FD[,-1]
traits.FD <- filter(traits.FD, row.names(traits.FD) %in% colnames(df.RA.wide))
df.RA.wide <- df.RA.wide[, colnames(df.RA.wide) %in% row.names(traits.FD)]
traits.FD <- traits.FD[order(row.names(traits.FD)),]

fd <- dbFD(as.matrix(traits.FD), as.matrix(df.RA.wide), stand.x = T, stand.FRic = T) 

saveRDS(fd, "Data/Long-Term-Datasets/FD-all.RDS")
# but hte functional traits for portal that I have only represent 73% of the community, wonder if this is changing things?
fd$FDiv
  # McLaughlin Jasper Ridge      Carrizo      Sonoran       Portal 
  #  0.9096691    0.8375451    0.8285393    0.7791593    0.8552983 

fd$FDis
 # McLaughlin Jasper Ridge      Carrizo      Sonoran       Portal 
 #    3.203720     2.991868     3.128179     2.657315     2.942931 

fd$RaoQ
  # McLaughlin Jasper Ridge      Carrizo      Sonoran       Portal 
  #  11.289721     9.145243    10.603856     7.711512     8.959819 

fd$FRic
#   McLaughlin Jasper Ridge      Carrizo      Sonoran       Portal 
# 1.783223e-01 1.336242e-05 1.749349e-02 9.683160e-04 5.374363e-06 
# still no idea what JR is so low!!

fd.all <- as.data.frame(fd$FRic)
fd.all$site <- row.names(fd.all)
colnames(fd.all)[1] <- "FRic"

fd.all <- merge(fd.all, meta[,c(1,19:22)], by.x = "site", by.y = "Dataset")

ggplot(fd.all, aes(x = log(AI.site), y = FRic)) +
  geom_point() +
  geom_smooth(method = "lm")

#### Time: FD ####
# functional diversity over time vs functional diversity over space? this could 
# tell us about spatial and temporal variability at each site
mcl <- readRDS("Data/Long-Term-Datasets/McLaughlin/McL_RA-time.RDS")
mcl$site <- "McLaughlin"
names(mcl)

jr <- readRDS("Data/Long-Term-Datasets/Jasper-Ridge/JR_RA-time.rds")
jr$site <- "Jasper Ridge"
colnames(jr)[2] <- "Year"

cp <- readRDS("Data/Long-Term-Datasets/Carrizo-Plain/Carrizo_RA-time.RDS")
cp$site <- "Carrizo"
names(cp)[1] <- "Code"
colnames(cp)[2] <- "Year"

sd <- readRDS("Data/Long-Term-Datasets/Sonoran-Desert/SD_RA-time.RDS")
sd$site <- "Sonoran"
names(sd)

portal <- readRDS("Data/Long-Term-Datasets/Portal/Portal_RA-time.RDS")
portal$site <- "Portal"
colnames(portal)[2] <- "Year"

df.RA <- rbind(mcl,jr,cp,sd,portal)

df.RA$Code <- recode_factor(df.RA$Code, 
                            FILCAL = "LOGFIL",
                            LYSARV = "ANAARV",
                            PHLGRA = "MICGRA",
                            PLAOVAI = "PLAINS")

df.RA.wide <- pivot_wider(sd[,-4], names_from = "Code", values_from = RA)
df.RA.wide <- as.data.frame(df.RA.wide)
row.names(df.RA.wide) <- df.RA.wide$Year
#df.RA.wide <- df.RA.wide[df.RA.wide$Year != 2006,] #portal, 1 species in 2006 but no traits for that species
df.RA.wide <- df.RA.wide[,-1]
df.RA.wide <- df.RA.wide[,order(colnames(df.RA.wide))]

traits.FD <- traits[,c(2,12:20,24,25,27)]
traits.FD <- as.data.frame(traits.FD)
row.names(traits.FD) <- traits.FD$new.code
traits.FD <- traits.FD[,-1]
traits.FD <- filter(traits.FD, row.names(traits.FD) %in% colnames(df.RA.wide))
df.RA.wide <- df.RA.wide[, colnames(df.RA.wide) %in% row.names(traits.FD)]
traits.FD <- traits.FD[order(row.names(traits.FD)),]
fd.sd <- dbFD(as.matrix(traits.FD), as.matrix(df.RA.wide), stand.FRic = T, stand.x = T)
fd.jr$FRic
mean(fd.cp$FRic)
fd.portal$FRic
fd.sd$FRic

fd.mcl.div <- as.data.frame(fd.mcl$FRic)
fd.mcl.div$Year <- row.names(fd.mcl.div)
fd.mcl.div$site <- "McLaughlin"
colnames(fd.mcl.div)[1] <- "FRic"

fd.portal.div <- as.data.frame(fd.portal$FRic)
fd.portal.div$Year <- row.names(fd.portal.div)
fd.portal.div$site <- "Portal"
colnames(fd.portal.div)[1] <- "FRic"

fd.cp.div <- as.data.frame(fd.cp$FRic)
fd.cp.div$Year <- row.names(fd.cp.div)
fd.cp.div$site <- "Carrizo"
colnames(fd.cp.div)[1] <- "FRic"

fd.sd.div <- as.data.frame(fd.sd$FRic)
fd.sd.div$Year <- row.names(fd.sd.div)
fd.sd.div$site <- "Sonoran"
colnames(fd.sd.div)[1] <- "FRic"

fd.jr.div <- as.data.frame(fd.jr$FRic)
fd.jr.div$Year <- row.names(fd.jr.div)
fd.jr.div$site <- "Jasper"
colnames(fd.jr.div)[1] <- "FRic"

fd <- rbind(fd.mcl.div, fd.portal.div, fd.cp.div, fd.sd.div, fd.jr.div)

ggplot(fd, aes(x = Year, y = FRic, col = site)) +
  geom_point() +
  facet_wrap(~site) # i dont get how these estimates of FRic are so wildly different from the site average FRic... by orders of magnitude

# the more the high abundances are greater than the mean, the higher FDiv is (villeger)


#### Space: FD ####

#### Time-Space: FD ####
# functional diversity over time vs functional diversity over space? this could 
# tell us about spatial and temporal variability at each site
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
fd.all.ts <- dbFD(as.matrix(traits.FD), as.matrix(df.RA.wide), w.abun = T, stand.x = T, stand.FRic = T) # WHY WERE SO MANY AXES REMOVED


## Functional richness ##
tmp <- as.data.frame(fd.all.ts$FRic)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "FRic"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp.sum <- tmp %>% 
  group_by(site, AI.site) %>%
  summarize(FRic.mean = mean(FRic, na.rm = T),
            FRic.se = calcSE(FRic))

ggplot(tmp.sum, aes(x = AI.site, y = FRic.mean)) +
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

## Rao's Q
tmp <- as.data.frame(fd.all.ts$RaoQ)
tmp$site <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[1]))
tmp$year <- sapply(strsplit(row.names(tmp), split='.', fixed=TRUE), function(x) (x[2]))
colnames(tmp)[1] <- "RaoQ"

tmp <- merge(tmp, meta[,c(1,22)], by.x = "site", by.y = "Dataset")

tmp.sum <- tmp %>% 
  group_by(site, AI.site) %>%
  summarize(RaoQ.mean = mean(RaoQ, na.rm = T),
            RaoQ.se = calcSE(RaoQ))

ggplot(tmp.sum, aes(x = AI.site, y = RaoQ.mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymax = RaoQ.mean + RaoQ.se, ymin = RaoQ.mean - RaoQ.se), width = 0.001)


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


# fd.mcl.div <- as.data.frame(fd.mcl$FRic)
# fd.mcl.div$Year <- row.names(fd.mcl.div)
# fd.mcl.div$site <- "McLaughlin"
# colnames(fd.mcl.div)[1] <- "FRic"
# 
# fd.portal.div <- as.data.frame(fd.portal$FRic)
# fd.portal.div$Year <- row.names(fd.portal.div)
# fd.portal.div$site <- "Portal"
# colnames(fd.portal.div)[1] <- "FRic"
# 
# fd.cp.div <- as.data.frame(fd.cp$FRic)
# fd.cp.div$Year <- row.names(fd.cp.div)
# fd.cp.div$site <- "Carrizo"
# colnames(fd.cp.div)[1] <- "FRic"
# 
# fd.sd.div <- as.data.frame(fd.sd$FRic)
# fd.sd.div$Year <- row.names(fd.sd.div)
# fd.sd.div$site <- "Sonoran"
# colnames(fd.sd.div)[1] <- "FRic"
# 
# fd.jr.div <- as.data.frame(fd.jr$FRic)
# fd.jr.div$Year <- row.names(fd.jr.div)
# fd.jr.div$site <- "Jasper"
# colnames(fd.jr.div)[1] <- "FRic"
# 
# fd <- rbind(fd.mcl.div, fd.portal.div, fd.cp.div, fd.sd.div, fd.jr.div)
# 
# ggplot(fd, aes(x = Year, y = FRic, col = site)) +
#   geom_point() +
#   facet_wrap(~site) # i dont get how these estimates of FRic are so wildly different from the site average FRic... by orders of magnitude
# 
# # the more the high abundances are greater than the mean, the higher FDiv is (villeger)