# HMM traits for CNPS
rm(list=ls())
library(lme4)
library(vegan)
library(Rmisc)
library(lmerTest)
library(nationalparkcolors)

#### Seed Traits ####
species <- read.csv("Data/20211001_Full-Species-List.csv")
traits <- read.csv("Data/20221011_Seed-Traits_cleaning.csv")
calipc <- read.csv("Collections/cal-ipc_inventory.csv")

# Normalize
traits$set.time.mpsec.sqrt <- sqrt(traits$set.time.mpsec)
traits$shape.asinsqrt <- asin(sqrt(traits$shape))
traits$mass.mg.log <- log(traits$mass.mg)
traits$height.cm.log <- log(traits$height.cm)
traits$wing.loading.log <- log(traits$wing.loading)
traits$cn <- ifelse(is.na(traits$cn), mean(traits$cn, na.rm = T), traits$cn)
traits$cn.sqrt <- sqrt(traits$cn)

traits$prop.C <- ifelse(is.na(traits$prop.C), mean(traits$prop.C, na.rm = T), traits$prop.C)
traits$C.prop.asinsqrt <- asin(sqrt(traits$prop.C))

traits$prop.N <- ifelse(is.na(traits$prop.N), mean(traits$prop.N, na.rm = T), traits$prop.N)
traits$N.prop.log <- log(traits$prop.N)

#traits$ldd <- hist(log(traits$ldd))
traits$size.log <- log10(traits$size)
traits$E.S <- ifelse(is.na(traits$E.S), mean(traits$E.S, na.rm = T), traits$E.S)
traits$E.S.asinsqrt <- asin(sqrt(traits$E.S))
traits$coat.perm <- ifelse(is.na(traits$coat.perm), mean(traits$coat.perm, na.rm = T), traits$coat.perm)
traits$coat.perm.log <- log(traits$coat.perm)

traits$coat.thick <- ifelse(is.na(traits$coat.thick), mean(traits$coat.thick, na.rm = T), traits$coat.thick)
traits$coat.thick.log <- log(traits$coat.thick)


traits <- merge(traits, calipc[,c(1,3)], by.x = "Species", by.y = "Scientific.name", all.x = T, all.y = F)

traits <- merge(traits, species[,c(3,7,8,9)], by = "Species")
traits$FunGroup <- paste(traits$nat.inv, traits$group, sep = " ")
traits[is.na(traits$Rating),]$Rating <- "None"
traits[traits$nat.inv == "native",]$Rating <- "Native"

disp <- c("Species", "code", "family", "group", "nat.inv", "FunGroup", "Rating", "appendage.type", "appendage", "ldd",  "set.time.mpsec.sqrt", "shape.asinsqrt", "mass.mg.log",  "height.cm.log", "wing.loading.log", "cn.sqrt", "size.log") #, "E.S.asinsqrt", "coat.perm.log", "coat.thick.log"

disp.pca <- c("ldd", "set.time.mpsec.sqrt", "shape.asinsqrt", "mass.mg.log",  "height.cm.log", "wing.loading.log", "cn.sqrt", "size.log") #, "E.S.asinsqrt", "coat.perm.log", "coat.thick.log"

no.disp <- c("Species", "code", "family", "group", "nat.inv", "FunGroup", "Rating", "appendage.type", "appendage", "ldd2", "set.time.mpsec.sqrt", "shape.asinsqrt", "mass.mg.log",  "height.cm.log", "wing.loading.log", "C.prop.asinsqrt", "N.prop.log",  "size.log") #, "E.S.asinsqrt", "coat.perm.log", "coat.thick.log"

no.disp.pca <- c("set.time.mpsec.sqrt", "shape.asinsqrt", "mass.mg.log",  "height.cm.log", "wing.loading.log",  "C.prop.asinsqrt", "N.prop.log",  "size.log") #, "E.S.asinsqrt", "coat.perm.log", "coat.thick.log"

#### .___PCA ALL ####

# with dispersal
trait.pca <- traits[complete.cases(traits$ldd), disp]

#trait.pca[trait.pca$Species == "Brassica tournefortii",]$ldd <- 3 ## BRASSICA TOUR HAS LDD, talk to Jutta, grey literature?
pca <- prcomp(trait.pca[, disp.pca], scale = T)
summary(pca)


trait.pca <- cbind(trait.pca, pca$x[,1:4])


autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "Rating", shape = "Rating") +
  theme_classic() +
  #geom_text(aes(label = code, col = Rating)) +
  #stat_ellipse(aes(group = group)) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "grey", "darkgoldenrod3"))
  #scale_color_manual(values = pal[c(1,3,4,5)])

autoplot(pca, x = 3, y = 4, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "Rating", shape = "Rating") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()) +
   #stat_ellipse(aes(group = FunGroup, col = FunGroup)) + 
   scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "grey", "darkgoldenrod3"))

# without dispersal
trait.pca <- traits[, no.disp]

pca <- prcomp(trait.pca[, no.disp.pca], scale = T)
summary(pca)

trait.pca <- cbind(trait.pca, pca$x[,1:4])

#pal <- park_palette("Redwoods")

no.data <- autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = F, label = F, col = "white") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank(),
    legend.position = "none"
  ) 

ggsave("/Users/Marina/Documents/USDA-PostDoc/Conferences/CNPS 2022/PCA_no_data.jpeg", no.data, units = "in", width = 3.5, height = 3, dpi = 600)

data.nat.inv <- autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = F, label = F, col = "nat.inv", shape = "nat.inv", size = 2) +
  #geom_point(aes(col = nat.inv, shape = nat.inv), size = 3) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  scale_color_manual(values = c("#1B9E77", "grey")) +
  scale_shape_manual(values = c(16, 3))

ggsave("/Users/Marina/Documents/USDA-PostDoc/Conferences/CNPS 2022/PCA_data_natinv.jpeg", data.nat.inv, units = "in", width = 3.5, height = 3, dpi = 600)

data <- autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = F, label = F, col = "Rating", shape = "Rating", size = 2) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "grey", "darkgoldenrod3"))
  
ggsave("/Users/Marina/Documents/USDA-PostDoc/Conferences/CNPS 2022/PCA_data.jpeg", data, units = "in", width = 3.5, height = 3, dpi = 600)

autoplot(pca, x = 3, y = 4, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "Rating", shape = "Rating") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "grey", "darkgoldenrod3"))

#### .___PCA Forbs ####
# disp
trait.pca.forb <- traits[traits$group == "forb", disp]
trait.pca.forb <- trait.pca.forb[complete.cases(trait.pca.forb),] 
pca.forb <- prcomp(trait.pca.forb[, disp.pca], scale = T)
summary(pca.forb)

trait.pca.forb <- cbind(trait.pca.forb, pca.forb$x[,1:4])


autoplot(pca.forb, x = 1, y = 2, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "Rating", shape = "Rating") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "grey", "darkgoldenrod3"))

autoplot(pca.forb, x = 3, y = 4, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "FunGroup") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
   ) #+
  # scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "grey", "darkgoldenrod3"))

# no disp
trait.pca.forb <- traits[traits$group == "forb", no.disp]
trait.pca.forb <- trait.pca.forb[complete.cases(trait.pca.forb),] 
pca.forb <- prcomp(trait.pca.forb[, no.disp.pca], scale = T)
summary(pca.forb)

trait.pca.forb <- cbind(trait.pca.forb, pca.forb$x[,1:4])


autoplot(pca.forb, x = 1, y = 2, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "Rating", shape = "Rating") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "grey", "darkgoldenrod3"))

autoplot(pca.forb, x = 3, y = 4, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "Rating", shape = "Rating") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "grey", "darkgoldenrod3"))

#### .___PCA McLaughlin only ####
mcl <- read.csv("/Users/Marina/Google Drive/02_McLaughlin_80Sites_Organized/Data-Storage-Files/Core_Community_Data2019.csv")

mcl.species <- unique(mcl$Species_Name)
mcl.species.old <- unique(mcl$Species_Name_J12)

traits.mcl <- traits[traits$Species %in% mcl.species,]
traits.mcl.old <- traits[traits$Species %in% mcl.species.old,]
traits.mcl <- unique(rbind(traits.mcl, traits.mcl.old))

# with dispersal
trait.pca.mcl <- traits.mcl[complete.cases(traits.mcl$ldd), disp]

#trait.pca[trait.pca$Species == "Brassica tournefortii",]$ldd <- 3 ## BRASSICA TOUR HAS LDD, talk to Jutta, grey literature?
pca <- prcomp(trait.pca.mcl[, disp.pca], scale = T)
summary(pca)


trait.pca.mcl <- cbind(trait.pca.mcl, pca$x[,1:4])

autoplot(pca, x = 1, y = 2, data = trait.pca.mcl, frame = F, loadings = T, loadings.label = T, label = F, col = "FunGroup", shape = "Rating") +
  theme_classic() +
  #geom_text(aes(label = code, col = Rating)) +
  #stat_ellipse(aes(group = group)) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "grey", "darkgoldenrod3"))

autoplot(pca, x = 3, y = 4, data = trait.pca.mcl, frame = F, loadings = T, loadings.label = T, label = F, col = "FunGroup", shape = "Rating") +
  theme_classic() +
  #geom_text(aes(label = code, col = Rating)) +
  #stat_ellipse(aes(group = group)) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "grey", "darkgoldenrod3"))

# without dispersal
trait.pca.mcl <- traits.mcl[, no.disp]

pca <- prcomp(trait.pca.mcl[, no.disp.pca], scale = T)
summary(pca)

trait.pca.mcl <- cbind(trait.pca.mcl, pca$x[,1:4])

autoplot(pca, x = 1, y = 2, data = trait.pca.mcl, frame = F, loadings = T, loadings.label = T, label = F, col = "FunGroup", shape = "Rating") +
  theme_classic() +
  #geom_text(aes(label = code, col = Rating)) +
  #stat_ellipse(aes(group = group)) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "grey", "darkgoldenrod3"))

autoplot(pca, x = 3, y = 4, data = trait.pca.mcl, frame = F, loadings = T, loadings.label = T, label = F, col = "FunGroup", shape = "Rating") +
  theme_classic() +
  #geom_text(aes(label = code, col = Rating)) +
  #stat_ellipse(aes(group = group)) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "grey", "darkgoldenrod3"))

#### HMMs ####
# mcl.df <- readRDS("HMMs/McLaughlin/mcl_HMM_quad_all_40.RDS")
# mcl.df$trt <- "all"
# 
# mcl.df2 <- readRDS("HMMs/McLaughlin/mcl_HMM_all-trt_quad_40.RDS")
# mcl.df <- rbind(mcl.df, mcl.df2)

#mcl.df <- readRDS("HMMs/McLaughlin/mcl_HMM_transect_10.RDS")
mcl.df <- readRDS("HMMs/McLaughlin/mcl_HMM_quad_all_40.RDS")

#### .___All Fungroups ####
#CNPS <- merge(mcl.df[mcl.df$trt == "all",], trait.pca[,-6], by.x = "Code", by.y = "code")

CNPS <- merge(mcl.df, trait.pca.mcl[,-6], by.x = "Code", by.y = "code")

CNPS <- filter(CNPS, iter < 150, Species_Name != "Phacelia parryi", Species_Name != "Trifolium microdon")

CNPS$FunGroup <- recode(CNPS$FunGroup, "Exotic Grass" = "Non-Native Grass", "Exotic Forb" = "Non-Native Forb")
  
test <- ggplot(CNPS[CNPS$FunGroup != "Native Grass",], aes(x = FunGroup, y = c, col = FunGroup)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_jitter(size = 2.5) + 
  #geom_text(aes(label = Code)) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
    scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(y = "Probability of Spatial Dispersal")
  
ggsave("/Users/Marina/Documents/USDA-PostDoc/Conferences/CNPS 2022/spatial.jpeg", test, units = "in", width = 6, height = 5, dpi = 600)

test2 <- ggplot(CNPS[CNPS$FunGroup != "Native Grass",], aes(x = FunGroup, y = s, col = FunGroup)) +
  geom_boxplot(outlier.shape = NA, size = 1) +
  geom_jitter(size = 2.5) + 
  #geom_text(aes(label = Code)) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(y = "Probability of Temporal Dispersal")

  
ggsave("/Users/Marina/Documents/USDA-PostDoc/Conferences/CNPS 2022/temporal.jpeg", test2, units = "in", width = 6, height = 5, dpi = 600)

ggplot(CNPS[CNPS$FunGroup != "Native Grass",], aes(x = FunGroup, y = g, col = FunGroup)) +
  geom_boxplot() +
  geom_jitter() + 
  #geom_text(aes(label = Code)) +
  theme_classic() 

CNPS.sum.s <- summarySE(CNPS[CNPS$FunGroup != "Native Grass",], measurevar = "s", groupvars = "FunGroup")
CNPS.sum.c <- summarySE(CNPS[CNPS$FunGroup != "Native Grass",], measurevar = "c", groupvars = "FunGroup")

ggplot(CNPS.sum.c, aes(x = FunGroup, y = c)) +
  geom_point() +
  geom_errorbar(aes(ymin = c - se, ymax = c + se), width = 0.1) 

ggplot(CNPS.sum.s, aes(x = FunGroup, y = s)) +
  geom_point() +
  geom_errorbar(aes(ymin = s - se, ymax = s + se), width = 0.1) 

ggplot(CNPS[CNPS$n.plots > 30,], aes(x = s, y = c, col = PC1)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  geom_text(aes(label = Code)) +
  theme_classic() +
  #scale_colour_gradient(low = "red", high = "blue") 
  scale_color_viridis_c() +
  labs(x = "Probability of Seed Survival", y = "Probability of Colonization")
  # scale_colour_gradient2(low = "red", mid = "gray",
  # high = "blue", midpoint = mean(CNPS$PC1))

ggplot(CNPS[CNPS$FunGroup != "Native Grass",], aes(x = s, y = c, col = FunGroup)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  #geom_text(aes(label = Code)) +
  #facet_wrap(~FunGroup) +
  theme_classic() +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3"))
  
  #scale_colour_gradient(low = "red", high = "blue")  
  #scale_colour_gradient2(low = "red", mid = "gray",
  #high = "blue", midpoint = mean(CNPS$PC2))

main <- ggplot(CNPS[CNPS$n.plots > 20 & CNPS$Code != "TRIMIC" & CNPS$Code != "EROBOT" ,], aes(x = s, y = c, col = PC2)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point(size = 2.5) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
  labs(x = "Probability of Seed Survival", y = "Probability of Colonization") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  scale_color_viridis_c()

ggsave("/Users/Marina/Documents/USDA-PostDoc/Conferences/CNPS 2022/tradeoff.jpeg", main, units = "in", width = 5, height = 4, dpi = 600)

m1 <- lm(s ~ PC2, data = CNPS[CNPS$n.plots > 20,])
summary(m1) # native forbs are shorter (shorter disp distances) and have higher seed coat permeability 

m2 <- lm(s ~ PC2, data = CNPS)
summary(m2)

m2 <- lm(c ~ s : PC2, data = CNPS)
summary(m2) # grasses have higher colonization and lower survival than native forbs


ggplot(CNPS, aes(x = s, y = c, col = PC4)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_point() +
  geom_text(aes(label = Code)) +
  #facet_wrap(~FunGroup) +
  theme_classic() +
  #scale_colour_gradient(low = "red", high = "blue")  
  scale_colour_gradient2(low = "red", mid = "gray",
  high = "blue", midpoint = mean(CNPS$PC4))

# primarily driven by exotic grasses
ggplot(CNPS, aes(x = s, y = c, col = PC2)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  geom_text(aes(label = Code)) +
  facet_wrap(~FunGroup) +
  theme_bw() +
  #scale_colour_gradient(low = "red", high = "blue")  
  scale_colour_gradient2(low = "red", mid = "gray",
  high = "blue", midpoint = mean(CNPS$PC2))

CNPS.sum <- CNPS %>% 
  dplyr::filter(FunGroup != "Native Grass") %>%
  dplyr::group_by(FunGroup) %>%
  dplyr::summarize(PC2.mean = mean(PC2), 
            PC2.se = sd(PC2)/sqrt(length((PC2))),
            PC3.mean = mean(PC3), 
            PC3.se = sd(PC3)/sqrt(length((PC3))),
            PC4.mean = mean(PC4), 
            PC4.se = sd(PC4)/sqrt(length((PC4))),
            PC1.mean = mean(PC1), 
            PC1.se = sd(PC1)/sqrt(length((PC1))),
            c.mean = mean(c),
            c.se = sd(c)/sqrt(length((c))),
            s.mean = mean(s),
            s.se = sd(s)/sqrt(length((s))))
  
ggplot(CNPS.sum, aes(x = PC2.mean, y = c.mean, col = FunGroup, group = FunGroup)) +
  geom_point() +
  geom_errorbar(aes(ymin = c.mean-c.se, ymax = c.mean+c.se)) + 
  geom_errorbarh(aes(xmin = PC2.mean - PC2.se, xmax =  PC2.mean + PC2.se)) + 
  theme_classic()

ggplot(CNPS.sum, aes(x = PC2.mean, y = s.mean, col = FunGroup, group = FunGroup)) +
  geom_point() +
  geom_errorbar(aes(ymin = s.mean - s.se, ymax = s.mean + s.se)) + 
  geom_errorbarh(aes(xmin = PC2.mean - PC2.se, xmax =  PC2.mean + PC2.se)) + 
  theme_classic()

ggplot(CNPS.sum, aes(x = PC4.mean, y = c.mean, col = FunGroup, group = FunGroup)) +
  geom_point() +
  geom_errorbar(aes(ymin = c.mean-c.se, ymax = c.mean+c.se)) + 
  geom_errorbarh(aes(xmin = PC4.mean - PC4.se, xmax =  PC4.mean + PC4.se)) + 
  theme_classic()

ggplot(CNPS.sum, aes(x = PC4.mean, y = s.mean, col = FunGroup, group = FunGroup)) +
  geom_point() +
  geom_errorbar(aes(ymin = s.mean - s.se, ymax = s.mean + s.se)) + 
  geom_errorbarh(aes(xmin = PC4.mean - PC4.se, xmax =  PC4.mean + PC4.se)) + 
  theme_classic()

#### .___Forbs only ####
CNPS.forb <- merge(mcl.df, trait.pca.forb[,-6], by.x = "Code", by.y = "code")

CNPS.forb <- filter(CNPS.forb, iter < 150, Species_Name != "Phacelia parryi", Species_Name != "Trifolium microdon", trt == "all")

ggplot(CNPS.forb, aes(x = s, y = c, col = PC1)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  geom_text(aes(label = Code)) +
  facet_wrap(~FunGroup) +
  theme_bw() +
  scale_colour_gradient2(low = "red", mid = "gray",
  high = "blue", midpoint = mean(CNPS.forb$PC1))

ggplot(CNPS.forb, aes(x = s, y = c, col = PC2)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  geom_text(aes(label = Code)) +
  facet_wrap(~FunGroup) +
  theme_bw() +
  scale_colour_gradient2(low = "red", mid = "gray",
  high = "blue", midpoint = mean(CNPS.forb$PC2))

ggplot(CNPS.forb, aes(x = s, y = c, col = PC3)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  geom_text(aes(label = Code)) +
  facet_wrap(~FunGroup) +
  theme_bw() +
  scale_colour_gradient2(low = "red", mid = "gray",
  high = "blue", midpoint = mean(CNPS.forb$PC3))

ggplot(CNPS.forb, aes(x = s, y = c, col = PC4)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  geom_text(aes(label = Code)) +
  facet_wrap(~FunGroup) +
  theme_bw() +
  scale_colour_gradient2(low = "red", mid = "gray",
  high = "blue", midpoint = mean(CNPS.forb$PC4))

ggplot(CNPS.forb, aes(x = PC2, y = c, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  #geom_text(aes(label = Code)) +
  theme_bw()

ggplot(CNPS.forb, aes(x = PC2, y = s, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  #geom_text(aes(label = Code)) +
  theme_bw()

ggplot(CNPS.forb, aes(x = s, y = c, col = PC2)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  #geom_text(aes(label = Code)) +
  facet_wrap(~FunGroup) +
  theme_bw() +
  scale_colour_gradient(low = "red", high = "blue")

#### MegaComp ####
mc.bio <- read.csv("Data/biomass-CNPS.csv")
mc.bio$total.biomass.g <- trimws(mc.bio$total.biomass.g)
mc.bio$bkgrd <- trimws(mc.bio$bkgrd)

mc.bio[mc.bio$phyto == "ACAM",]$phyto <- "ACMAME"
mc.bio[mc.bio$phyto == "ANAR",]$phyto <- "ANAARV"
mc.bio[mc.bio$phyto == "AVBA",]$phyto <- "AVEBAR"
mc.bio[mc.bio$phyto == "BRHO",]$phyto <- "BROHOR"
mc.bio[mc.bio$phyto == "GITR",]$phyto <- "GILTRI"
mc.bio[mc.bio$phyto == "MICA",]$phyto <- "MICCAL"
mc.bio[mc.bio$phyto == "THIR-I",]$phyto <- "TRIHIR"
mc.bio[mc.bio$phyto == "THIR-I ",]$phyto <- "TRIHIR"
mc.bio[mc.bio$phyto == "TWIL-I",]$phyto <- "TRIWIL"
mc.bio[mc.bio$phyto == "LENI",]$phyto <- "LEPNIT"
mc.bio[mc.bio$phyto == "PLER",]$phyto <- "PLAERE"

mc.bio[mc.bio$bkgrd == "ACAM",]$bkgrd <- "ACMAME"
mc.bio[mc.bio$bkgrd == "ANAR",]$bkgrd <- "ANAARV"
mc.bio[mc.bio$bkgrd == "AVBA",]$bkgrd <- "AVEBAR"
mc.bio[mc.bio$bkgrd == "BRHO",]$bkgrd <- "BROHOR"
mc.bio[mc.bio$bkgrd == "GITR",]$bkgrd <- "GILTRI"
mc.bio[mc.bio$bkgrd == "MICA",]$bkgrd <- "MICCAL"
mc.bio[mc.bio$bkgrd == "THIR-I",]$bkgrd <- "TRIHIR"
mc.bio[mc.bio$bkgrd == "TWIL-I",]$bkgrd <- "TRIWIL"
mc.bio[mc.bio$bkgrd == "LENI",]$bkgrd <- "LEPNIT"
mc.bio[mc.bio$bkgrd == "PLER",]$bkgrd <- "PLAERE"
mc.bio[mc.bio$bkgrd == "AMME",]$bkgrd <- "AMSMEN"
mc.bio[mc.bio$bkgrd == "BRNI",]$bkgrd <- "BRANIG"
mc.bio[mc.bio$bkgrd == "CESO",]$bkgrd <- "CENSOL"
mc.bio[mc.bio$bkgrd == "CLPU",]$bkgrd <- "CLAPUR"
mc.bio[mc.bio$bkgrd == "LOMU",]$bkgrd <- "FESPER"
mc.bio[mc.bio$bkgrd == "MAEL",]$bkgrd <- "MADELE"
mc.bio[mc.bio$bkgrd == "PLNO",]$bkgrd <- "PLANOT"
mc.bio[mc.bio$bkgrd == "TACA",]$bkgrd <- "ELYCAP"

mc.bio$bkgrd.group <- ifelse(mc.bio$bkgrd %in% species[species$nat.inv == "native",]$code, "native", "invasive")
mc.bio$bkgrd.group <- ifelse(mc.bio$bkgrd == "Control", "None", mc.bio$bkgrd.group)

mc.bio <- merge(mc.bio, trait.pca, by.x = "phyto", by.y = "code", all.y = F)

treat <- data.frame(block = unique(mc.bio$block), trt = c("D", "D", "D", "C", "D", "C", "C", "D", "D", "C", "C"))
mc.bio <- merge(mc.bio, treat, by = "block")
str(mc.bio)

mc.bio$total.biomass.g <- as.numeric(mc.bio$total.biomass.g)/mc.bio$phyto.n.indiv

mc.bio$seed.num <- mc.bio$seed.num/mc.bio$phyto.n.indiv

mc.bio <- filter(mc.bio, !bkgrd %in% c("THIR-U", "TWIL-U", "VIVI", "ERBO"))

ggplot(mc.bio, aes(x = PC4, y = log(seed.num), col = trt, group = trt)) +
  geom_point() +
  facet_wrap(~dens + bkgrd.group) + 
  geom_smooth(method = "lm", se = F)

m.bio <- lmer(log(seed.num) ~ PC2 + (1|block), data = mc.bio[mc.bio$seed.num > 0 & mc.bio$nat.inv == "native",])
plot(fitted(m.bio), resid(m.bio))
qqnorm(resid(m.bio))
qqline(resid(m.bio))
summary(m.bio)

mc.bio.sum <- summarySE(mc.bio, measurevar = "seed.num", groupvars = c("nat.inv", "bkgrd.group", "PC2", "PC3", "PC4"), na.rm = T)

ggplot(mc.bio[mc.bio$bkgrd.group == "native" & mc.bio$nat.inv == "invasive",], aes(x = PC2, y = log(total.biomass.g))) +
  geom_point() + 
  #geom_line() + 
  #geom_errorbar(aes(ymin = log(seed.num - se), ymax = log(seed.num + se)), width = 0.1)

ggplot(mc.bio.sum[mc.bio.sum$bkgrd.group == "None",], aes(x = PC3, y = log(seed.num))) +
  geom_point() 
  #geom_line() + 
  #geom_errorbar(aes(ymin = log(seed.num - se), ymax = log(seed.num + se)), width = 0.1)

ggplot(mc.bio.sum[mc.bio.sum$bkgrd.group == "None",], aes(x = PC4, y = log(seed.num))) +
  geom_point() + 
  #geom_line() + 
  geom_errorbar(aes(ymin = log(seed.num - se), ymax = log(seed.num + se)), width = 0.1)

mc.bio.sum <- summarySE(mc.bio, measurevar = "seed.num", groupvars = c("PC2", "bkgrd.group"), na.rm = T)

ggplot(mc.bio.sum, aes(x = PC2, y = log(seed.num), col = bkgrd.group, group = bkgrd.group)) +
  geom_point() +
  geom_errorbar(aes(ymin = log(seed.num - se), ymax = log(seed.num + se)), width = 0.1)  

traits.mc <- read.csv("Data/20210610_Seeds_All-Accessions.csv")
traits.mc <- filter(traits.mc, mega.comp == "x")



test <- filter(trait.pca, Species %in% traits.mc$Species)

#### Germination ####
germ <- read.csv("Data/20220218_Germination-Data_full.csv")
viab <- read.csv("Data/CURE_Viability_entered.csv")
viab$p.viab <- viab$Viable/viab$Seeds.left

germ <- filter(germ, Species != "EROBOT")
test <- lmer(p.germ ~ Group * Treatment + (1|Species), data = germ)
plot(fitted(test), resid(test))
qqnorm(resid(test))
qqline(resid(test))
summary(test)
anova(test) #, group, treatment, and group*treatment effect

germ.sum <- summarySE(germ, measurevar = "p.germ", groupvars = c("Treatment", "Group"))

ggplot(germ.sum, aes(x = Treatment, y = p.germ, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = p.germ - se, ymax = p.germ + se, width = 0.5)) +
  facet_wrap(~Group)

germ <- merge(germ, trait.pca, by.x = "Species", by.y = "code", all.y = F)

ggplot(germ, aes(x = PC2, y = p.germ)) +
  geom_point() +
  facet_wrap(~Treatment)

germ.sum <- summarySE(germ, measurevar = "p.germ", groupvars = c("Treatment", "PC2", "Group", "Species"))

ggplot(germ.sum, aes(x = PC2, y = p.germ, col = Group)) +
  geom_point() + 
  geom_text(aes(label = Species)) +
  #geom_smooth(method = "lm", se = F) +
  facet_wrap(~Treatment)

test <- lmer(p.germ ~ PC2 * Treatment + (1|Species), data = germ)
plot(fitted(test), resid(test))
qqnorm(resid(test))
qqline(resid(test))
summary(test)
anova(test) #, group, treatment, and group*treatment effect

germ.sum2 <- summarySE(germ, measurevar = "p.germ", groupvars = c("Treatment", "shape.asinsqrt", "Group", "Species"))

ggplot(germ.sum2, aes(x = shape.asinsqrt, y = p.germ, col = Group)) +
  geom_point() +
  #geom_text(aes(label = Species)) +
  geom_errorbar(aes(ymin = p.germ - se, ymax = p.germ + se, width = 0.01)) +
  facet_wrap(~Treatment)

# could also look at seed traits WITH reproductive output averages since that is typically an important trait