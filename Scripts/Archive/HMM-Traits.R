#### HMM System Comparison ####
rm(list=ls())

calcSE<-function(x){
  x2<-na.omit(x)
  sd(x2)/sqrt(length(x2))
}


library(ggfortify)
library(Rmisc)
library(tidyverse)
library(gridExtra)
library(sjstats)
# library(devtools)
# devtools::install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)
library(ape)
library(picante)
library(phytools)
library(Rphylopars)

#### Prep Data ####
species <- read.csv("Data/20211001_Full-Species-List.csv")

# Carrizo #
#CP.df2 <- readRDS("Scripts/HMMs/Carrizo-Plain/CP_HMM_all_40.RDS")
CP.df <- readRDS("Scripts/HMMs/Carrizo-Plain/CP_HMM_gkr_40.RDS")
CP.df <- filter(CP.df, Code != "EROCIC")
CP.df$trt <- "all"

# Jasper Ridge #
jr.df <- readRDS("Scripts/HMMs/Jasper-Ridge/jr_HMM_05m_all.RDS") 
#jr.df5[jr.df5$trt == "grc",]$trt <- "all"
#jr.df <- jr1.grc[jr1.grc$trt == "all",]
jr.df$trt <- "all"
#jr.df1 <- readRDS("Scripts/HMMs/Jasper-Ridge/jr_HMM_1m.RDS")
#jr.df1[jr.df1$trt == "grc",]$trt <- "all"

# Sonoran Desert #
sd.df <- readRDS("Scripts/HMMs/Sonoran-Desert/SD_HMM_all-2001_01m2-only.RDS")
#sd.df <- readRDS("Scripts/HMMs/Sonoran-Desert/SD_HMM_all-1995.RDS")
sd.df$trt <- "all"
#sd.df <- filter(sd.df, n.plots >= 8)

# McLaughlin #
mcl.df <- readRDS("Scripts/HMMs/McLaughlin/20221212_mcl_HMM_quadrat.RDS")[,-14]
mcl.df$trt <- "all"
mcl.df <- filter(mcl.df, n.plots >= 80)

# Portal #
port.df <- readRDS("Scripts/HMMs/Portal/Port_HMM_all-combined_40.RDS")
port.df <- merge(port.df, species[,c(3,6)], by.x = "Code", by.y = "code")
port.df <- port.df[,c(16,1,4,5:15)]
colnames(port.df)[1] <- "Species_Name"

# Add sites #
sd.df$site <- "Sonoran Desert"
port.df$site <- "Chihuahuan Desert"
mcl.df$site <- "McLaughlin"
CP.df$site <- "Carrizo"
jr.df1$site <- "Jasper Ridge"
jr.df5$site <- "Jasper Ridge"
jr.df$site <- "Jasper Ridge"

# Add AI
sd.df$AI <- 0.1157
port.df$AI <- 0.1320
mcl.df$AI <- 0.5511
CP.df$AI <- 0.1883
jr.df1$AI <- 0.4011
jr.df5$AI <- 0.4011
jr.df$AI <- 0.4011

# bind together #
full <- rbind(CP.df, jr.df, mcl.df, sd.df, port.df)

#full <- rbind(CP.df, jr.df1[,-14], mcl.df, sd.df, port.df)

#rm(CP.df, jr.df, mcl.df, sd.df, port.df)


#### Seed Traits ####
traits <- read.csv("Data/20230607_Seed-Traits_clean_site.csv")
species.AI <- read.csv("Data/20230530_Seeds_All-Accessions.csv")

# AMSMEN has a suspiciously low prop.C compared to other Amsinckia species, replace with mean
traits[traits$Species == "Amsinckia menziesii",]$prop.C <- mean(traits$prop.C, na.rm = T)

traits <- filter(traits, new.code != "CUCPAL") # remove because ??
traits[is.na(traits$prop.C),]$prop.C <- mean(traits$prop.C, na.rm = T)
traits[is.na(traits$prop.N),]$prop.N <- mean(traits$prop.N, na.rm = T)
traits[is.na(traits$coat.perm.perc),]$coat.perm.perc <- mean(traits$coat.perm.perc, na.rm = T)

traits[is.na(traits$wing.loading),]$wing.loading <- mean(traits$wing.loading, na.rm = T)

# Normalize
hist(log(traits$wing.loading))
traits$wing.loading <- log(traits$wing.loading)

hist(sqrt(traits$prop.C))
traits$prop.C <- sqrt(traits$prop.C)

hist(log(traits$prop.N))
traits$prop.N <- log(traits$prop.N)

hist(log(traits$coat.perm.perc))
traits$coat.perm.perc <- log(traits$coat.perm.perc)

hist(log(traits$morph.mass.mg))
traits$morph.mass.mg <- log(traits$morph.mass.mg)

hist(log(traits$chem.mass.mg))
traits$chem.mass.mg <- log(traits$chem.mass.mg)

hist(log(traits$size.mm))
traits$size.mm <- log(traits$size.mm)

hist(traits$set.time.mpsec)

hist(sqrt(traits$shape))
traits$shape <- sqrt(traits$shape)

hist(log(traits$E.S))
traits$E.S <- log(traits$E.S)

traits$coat.thick.per.width <- traits$both.thick/traits$width

hist(log(traits$coat.thick.per.width))
traits$coat.thick.per.width <- log(traits$coat.thick.per.width)

traits$FunGroup <- paste(traits$nat.inv, traits$group, sep = " ")

traits.ID <- traits

traits <- traits %>%
  group_by(Species, new.code, family, group, nat.inv, FunGroup, appendage.type, appendage, ldd2)  %>%
  dplyr::summarize(across(morph.mass.mg:coat.thick.per.width, ~ mean(.x, na.rm = TRUE)))


trait.col <- c("Species", "new.code", "family", "group", "nat.inv", "FunGroup", "chem.mass.mg",  "shape", "size.mm", "appendage.type", "appendage", "prop.C", "prop.N", "set.time.mpsec", "height.cm", "ldd2", "coat.perm.perc", "E.S", "both.thick", "mucilage", "fleshy.end", "starchy.end", "coat.thick.per.width", "ldd", "AI", "AI.mean", "AI.IQR", "AI.25", "AI.75")

trait.pca.col <- c("chem.mass.mg",  "shape", 
"prop.C", "prop.N", "set.time.mpsec", "height.cm", "coat.perm.perc", "E.S", "size.mm", "coat.thick.per.width")

#### .___PCA ALL ####

trait.pca <- traits[, trait.col]

pca <- prcomp(trait.pca[, trait.pca.col], scale = T)
summary(pca)

trait.pca <- cbind(trait.pca, pca$x[,1:4])

autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "AI.mean") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_gradient2(low = "red", high = "blue", mid = "grey", midpoint = median(trait.pca$AI.mean))

autoplot(pca, x = 3, y = 4, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "ldd2") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  )

autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "system") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  )

autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "AI.mean") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_gradient2(low = "red", high = "blue", mid = "grey", midpoint = 0.3)


autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "ldd2") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  )

autoplot(pca, x = 3, y = 4, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "ldd2") +
  theme_classic() +
  #geom_text(aes(label = code, col = Rating)) +
  #stat_ellipse(aes(group = group)) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  )

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
trait.pca.forb <- traits[traits$group == "forb", trait.col]

pca.forb <- prcomp(trait.pca.forb[, trait.pca.col], scale = T)
summary(pca.forb)

trait.pca.forb <- cbind(trait.pca.forb, pca.forb$x[,1:4])


autoplot(pca.forb, x = 1, y = 2, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "ldd2") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  )

autoplot(pca.forb, x = 2, y = 3, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "ldd2") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  )

autoplot(pca.forb, x = 1, y = 2, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "FunGroup") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  )

autoplot(pca.forb, x = 3, y = 4, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "ldd2") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
   ) #+
  # scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "grey", "darkgoldenrod3"))

# hmm.trait <- merge(full, trait.pca, by.x = "Code", by.y = "code")
# 
# hmm.trait <- filter(hmm.trait, iter < 150, Species_Name != "Phacelia parryi", Species_Name != "Trifolium microdon") # need to double check that there aren't any cases like this, PHAPAR in the seed trait dataset matches with phalaris paradoxa in the  dataset, also microdon pulls up microcephalum

#### Seed traits vs aridity ####



#### Merge HMM and Seed traits ####
trait.pca$FunGroup <- ifelse(trait.pca$FunGroup == "native forb", "Native Forb", ifelse(trait.pca$FunGroup == "native grass", "Native Grass", trait.pca$FunGroup))

trait.pca$FunGroup <- ifelse(trait.pca$FunGroup == "invasive forb", "Exotic Forb", ifelse(trait.pca$FunGroup == "invasive grass", "Exotic Grass", trait.pca$FunGroup))

trait.pca.forb$FunGroup <- ifelse(trait.pca.forb$FunGroup == "native forb", "Native Forb", ifelse(trait.pca.forb$FunGroup == "native grass", "Native Grass", trait.pca.forb$FunGroup))

trait.pca.forb$FunGroup <- ifelse(trait.pca.forb$FunGroup == "invasive forb", "Exotic Forb", ifelse(trait.pca.forb$FunGroup == "invasive grass", "Exotic Grass", trait.pca.forb$FunGroup))

full <- filter(full, trt == "all")
hmm.trait <- merge(full, trait.pca.forb, by.x = c("Code", "FunGroup"), by.y = c("code", "FunGroup"))

hmm.trait <- filter(hmm.trait, iter < 150, Species_Name != "Phacelia parryi", Species_Name != "Trifolium microdon") 

# need to go back and take these species out of MCL see susan email about ID confusion
hmm.trait <- filter(hmm.trait, Code != "GERDIS", Code != "GERMOLL", Code != "SILGAL")

#### Graphs #### 

####.___PC1 ####
hmm.trait$system <- ifelse(hmm.trait$AI < 0.2, "desert", "grassland")

ggplot(hmm.trait[hmm.trait$iter < 150,], aes(x = s, y = c, col = system, group = system)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

ggplot(hmm.trait[hmm.trait$iter < 150 & hmm.trait$trt == "all",], aes(x = PC1, y = s, col = system, group = system)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

ggplot(hmm.trait[hmm.trait$iter < 150 & hmm.trait$trt == "all",], aes(x = PC1, y = c, col = system, group = system)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw()

####.___PC2 ####

ggplot(hmm.trait[hmm.trait$iter < 150 & hmm.trait$trt == "all",], aes(x = PC2, y = s, col = system, group = system)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw()

ggplot(hmm.trait[hmm.trait$iter < 150 & hmm.trait$trt == "all",], aes(x = PC2, y = c, col = system, group = system)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw()

####.___PC3####
ggplot(hmm.trait[hmm.trait$iter < 150 & hmm.trait$trt == "all",], aes(x = PC3, y = s, col = system, group = system)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw()

ggplot(hmm.trait[hmm.trait$iter < 150 & hmm.trait$trt == "all",], aes(x = PC3, y = c, col = system, group = system)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw()

#### Traits by system ####
ggplot(hmm.trait[hmm.trait$iter < 150,], aes(x = system, y = shape)) +
  geom_boxplot() 

ggplot(hmm.trait[hmm.trait$iter < 150,], aes(x = system, y = size.mm)) +
  geom_boxplot() 

ggplot(hmm.trait[hmm.trait$iter < 150,], aes(x = system, y = prop.C)) +
  geom_boxplot() 

ggplot(hmm.trait[hmm.trait$iter < 150,], aes(x = system, y = coat.thick.per.width)) +
  geom_boxplot() 

ggplot(traits[traits$group == "forb",], aes(x = system, y = shape)) +
  geom_boxplot() 

ggplot(traits[traits$group == "forb",], aes(x = system, y = size.mm)) +
  geom_boxplot() 

ggplot(traits[traits$group == "forb",], aes(x = system, y = prop.C)) +
  geom_boxplot() 

ggplot(traits[traits$group == "forb",], aes(x = system, y = coat.thick.per.width)) +
  geom_boxplot() 

#### Hier Part ####
survival <- c("chem.mass.mg",  "shape", "size.mm", "prop.C", "prop.N", "coat.perm.perc", "coat.thick.per.width")

colonization <- c("shape", "chem.mass.mg", "set.time.mpsec", "height.cm", "size.mm", "ldd2")

# pick traits without correlations but RDA might be able to handle them together

# how is this different from a hierarchical model? order matters
# conceptually it's fairly similar 

##### All Forbs ####

###### S: All Forbs ####

HPA.s <- rdacca.hp(
  dv = hmm.trait[hmm.trait$system == "grassland",]$s,
  iv = hmm.trait[hmm.trait$system == "grassland", colnames(hmm.trait) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation

plot.rdaccahp(HPA.s) 

HPA.s <- rdacca.hp(
  dv = hmm.trait[hmm.trait$system == "desert",]$s,
  iv = hmm.trait[hmm.trait$system == "desert", colnames(hmm.trait) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation

plot.rdaccahp(HPA.s) 

###### C: All Forbs ####

HPA.c <- rdacca.hp(
  dv = hmm.trait[hmm.trait$system == "grassland",]$c,
  iv = hmm.trait[hmm.trait$system == "grassland", colnames(hmm.trait) %in% colonization],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.c$Hier.part
#HPA.c$Var.part
HPA.c$Total_explained_variation

plot.rdaccahp(HPA.c) # size and shape

HPA.c <- rdacca.hp(
  dv = hmm.trait[hmm.trait$system == "desert",]$c,
  iv = hmm.trait[hmm.trait$system == "desert", colnames(hmm.trait) %in% colonization],
  method = "RDA",
  type = "adjR2",
  var.part = T
)

HPA.c$Hier.part
#HPA.c$Var.part
HPA.c$Total_explained_variation

plot.rdaccahp(HPA.c) # size and shape

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
traits.forbs <- traits[traits$group == "forb",]

# positive: mass, settling time, wing loading, coat thickness, 
# negative: prop.C

pairs(traits.forbs[, c(34,35,37, 10:15)], 
      lower.panel = panel.lmline, 
      upper.panel = panel.cor)

pairs(traits.forbs[, c(34,35,37, 16:20)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

pairs(traits.forbs[, c(34,35,37, 21:25)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

pairs(traits.forbs[, c(34,35,37, 26:30,39)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

traits$AI.group <- ifelse(traits$AI.mean <= 0.2, "arid", 
                          ifelse(traits$AI.mean > 0.2 & traits$AI.mean <= 0.5, "semi-arid", "dry sub-humid"))

traits$AI.group <- ifelse(traits$AI.mean > 0.65, "humid", traits$AI.group)

traits$AI.lowrange <- ifelse(traits$AI.25 <= 0.2, "arid", 
                          ifelse(traits$AI.25 > 0.2 & traits$AI.25 <= 0.5, "semi-arid", "dry sub-humid"))

traits$AI.lowrange <- ifelse(traits$AI.25 > 0.65, "humid", traits$AI.lowrange)

traits$AI.uprange <- ifelse(traits$AI.75 <= 0.2, "arid", 
                          ifelse(traits$AI.75 > 0.2 & traits$AI.75 <= 0.5, "semi-arid", "dry sub-humid"))

traits$AI.uprange <- ifelse(traits$AI.75 > 0.65, "humid", traits$AI.uprange)

write.csv(traits, "Data/test.csv", row.names = F)

# traits$AI.sum <- ifelse(traits$AI.group == "arid" & 
#                           traits$AI.uprange == "arid" &
#                           traits$AI.lowrange == "arid",
#                         "arid",
#                         )



# ok but shouldnt community weighted traits matter more? just because they are the same doesnt mean shit
''
traits <- read.csv("Data/test.csv")

traits$AI.sum <- factor(traits$AI.sum, levels = c("arid", "arid-preferred", "semi-arid preferred", "semi-arid"))

traits$AI.bi <- ifelse(traits$AI.sum == "arid" | traits$AI.sum == "arid-preferred", "arid", "semi-arid")

traits.sum <- traits %>%
  filter(group == "forb") %>%
  group_by(nat.inv, AI.bi) %>%
  dplyr::summarize(across(morph.mass.mg:coat.thick.per.width,
                          list(mean = ~mean(.x, na.rm = TRUE),
                               se = ~ calcSE(.x))))


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

ggplot(traits.sum, aes(x = AI.bi, y = shape_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = shape_mean - shape_se, ymax = shape_mean + shape_se, width = 0.1)) +
  facet_wrap(~nat.inv)

ggplot(traits.sum, aes(x = AI.bi, y = size.mm_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = size.mm_mean - size.mm_se, ymax = size.mm_mean + size.mm_se, width = 0.1)) +
  facet_wrap(~nat.inv)

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