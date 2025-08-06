### HMM Validation ###

rm(list=ls())

library(tidyverse)
library(Metrics)
library(ggplot2)
library(ggpubr)
library(ggfortify)

#### Sonoran Desert ####
##### Boots #####
boot <- readRDS("Scripts/HMMs/Sonoran-Desert/SD-species-boot.RDS")
trait <- read.csv("Data/Long-Term-Datasets/Sonoran-Desert/Species_list_2016_2.csv")
sd.meta <- read.csv("Data/Long-Term-Datasets/Sonoran-Desert/HMM-meta-SD.csv")

boot <- as.data.frame(boot)

boot <- boot %>% 
  dplyr::mutate(across(p0.boot:perc.con, ~as.numeric(.x)))

colnames(boot)[1] <- "code.old"
boot <- merge(boot, trait[,c(3,4)], by.x = "code.old", by.y = "Species.code")

boot  <- boot[,c(15,1:14)]

boot <- boot[boot$s.boot.sd < 0.1 & boot$c.boot.sd < 0.1 & boot$perc.con > 0.5,]

boot <- merge(boot, SD.df[,c(2,10:13)], by = "Code")

ggplot(boot[boot$max.years >= 3 & boot$n.plots >= 15,], aes(x = s.boot, y = c.boot)) +
  geom_smooth(method = "lm") +
  geom_text(aes(label = Code)) +
  geom_errorbar(aes(xmin = s.boot - s.boot.sd, xmax = s.boot + s.boot.sd)) +
  geom_errorbar(aes(ymin = c.boot - c.boot.sd, ymax = c.boot + c.boot.sd))


##### Patches #####
SD.sim.sp15 <- readRDS("Scripts/HMMs/Sonoran-Desert/SD-sim-species-validation-15.RDS")
SD.sim.sp20 <- readRDS("Scripts/HMMs/Sonoran-Desert/SD-sim-species-validation-20.RDS")
SD.sim.sp30 <- readRDS("Scripts/HMMs/Sonoran-Desert/SD-sim-species-validation-30.RDS")

SD.sim.sp <- rbind(SD.sim.sp15, SD.sim.sp20, SD.sim.sp30)

SD.sim.sp <- as.data.frame(SD.sim.sp)

SD.sim.sp <- SD.sim.sp %>% 
  dplyr::mutate(across(sim:iter, ~as.numeric(.x)))

colnames(SD.sim.sp)[4:9] <- paste(colnames(SD.sim.sp)[4:9], "est", sep=".")

SD.sim.sp <- merge(trait[, c(3,4)], SD.sim.sp, by.x = "Species.code", by.y = "Species_Name", all.x = F, all.y = T)

SD.sim.sp <- merge(SD.sim.sp, boot[,1:7], by = "Code", all.x = T, all.y = F)

SD.sim.sp <- merge(SD.sim.sp, SD.df[,c(2,10:13)], by = "Code")

SD.rmse <- SD.sim.sp %>%
  dplyr::group_by(sim, patches) %>%
  dplyr::filter(mean.years >= 2) %>%
  dplyr::summarize(s.bias = mean(abs(s.boot - mean(s.est))),
            s.rmse = rmse(s.boot, s.est),
            s.var = var(s.boot, s.est),
            c.bias = mean(abs(c.boot - mean(c.est))),
            c.rmse = rmse(c.boot, c.est),
            c.var = var(c.boot, g.est),
            g.bias = mean(abs(g.boot - mean(g.est))),
            g.rmse = rmse(g.boot, g.est),
            g.var = var(g.boot, g.est)
            )

a <- ggplot(SD.rmse, aes(x = s.rmse, group = patches, col = patches)) +
  geom_density(adjust = 2) +
  geom_vline(xintercept = mean(SD.rmse[SD.rmse$patches == 15,]$s.rmse), col = '#1b9e77', linetype = 2) +
  geom_vline(xintercept = mean(SD.rmse[SD.rmse$patches == 20,]$s.rmse), col = '#d95f02', linetype = 2) +
  geom_vline(xintercept = mean(SD.rmse[SD.rmse$patches == 30,]$s.rmse), col = '#7570b3', linetype = 2) +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) + theme_classic() +
  theme(
      legend.position = c(0.2,0.8)
  ) +
  labs(x = "RMSE of temporal dispersal estimates") +
  xlim(0,0.4)

b <- ggplot(SD.rmse, aes(x = c.rmse, group = patches, col = patches)) +
  geom_density(adjust = 1.5) +
  geom_vline(xintercept = mean(SD.rmse[SD.rmse$patches == 15,]$c.rmse), col = '#1b9e77', linetype = 2) +
  geom_vline(xintercept = mean(SD.rmse[SD.rmse$patches == 20,]$c.rmse), col = '#d95f02', linetype = 2) +
  geom_vline(xintercept = mean(SD.rmse[SD.rmse$patches == 30,]$c.rmse), col = '#7570b3', linetype = 2) +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) + theme_classic() +
  theme(
      legend.position = c(0.8,0.8)
  ) +
  labs(x = "RMSE of spatial dispersal estimates") +
  xlim(0,0.3)

SD.sims.fig <- ggarrange(a,b)

SD.sims.fig <- annotate_figure(SD.sims.fig, top = text_grob("Sonoran Desert (max 44 patches, 16 years)", color = "black", face = "bold", size = 16))

ggplot(SD.rmse, aes(x = s.bias, group = patches, col = patches)) +
  geom_density(adjust = 1.5)

ggplot(SD.rmse, aes(x = s.var, group = patches, col = patches)) +
  geom_density(adjust = 1.5)

#### Portal ####

##### Boots ####
boot.w <- readRDS("Scripts/HMMs/Portal/port-species-boot-winter.RDS")
boot.s <- readRDS("Scripts/HMMs/Portal/port-species-boot-summer.RDS")
boot.sw <- readRDS("Scripts/HMMs/Portal/port-species-boot-summer-winter.RDS")
trait.repo <- read.csv("Data/Long-Term-Datasets/Portal/repo-traits.csv")
trait.bio <- read.csv("Data/Long-Term-Datasets/Portal/biotime-portal-species.csv")
port.meta.w <- read.csv("Data/Long-Term-Datasets/Portal/HMM-meta-port-winter.csv")
port.meta.s <- read.csv("Data/Long-Term-Datasets/Portal/HMM-meta-port-summer.csv")
port.meta.sw <- read.csv("Data/Long-Term-Datasets/Portal/HMM-meta-port-summer-winter.csv")

traits <- merge(trait.repo, trait.bio, by = "Code", all = T)

boot.w <- as.data.frame(boot.w)
boot.s <- as.data.frame(boot.s)
boot.sw <- as.data.frame(boot.sw)
boot <- rbind(boot.w, boot.s, boot.sw)

boot <- boot %>% 
  dplyr::mutate(across(p0.boot:perc.con, ~as.numeric(.x)))

traits <- merge(trait.repo, trait.bio, by = "Code", all = T)

boot <- merge(traits[, c(1,2)], boot, by.x = "species", by.y = "Species_Name", all.x = F, all.y = T)

boot <- boot[,c(2:15)]

boot <- boot[boot$s.boot.sd < 0.1 & boot$c.boot.sd < 0.1 & boot$perc.con > 0.5,]

port.meta <- rbind(port.meta.w, port.meta.s, port.meta.sw)

boot <- merge(boot, port.meta[,2:8], by = "Code")

ggplot(boot, aes(x = s.boot, y = c.boot)) +
  geom_smooth(method = "lm") +
  geom_text(aes(label = Code)) +
  geom_errorbar(aes(xmin = s.boot - s.boot.sd, xmax = s.boot + s.boot.sd)) +
  geom_errorbar(aes(ymin = c.boot - c.boot.sd, ymax = c.boot + c.boot.sd))


##### Patches #####
port.sim.sp15 <- readRDS("Scripts/HMMs/Portal/port-sim-species-validation-15.RDS")
port.sim.sp20 <- readRDS("Scripts/HMMs/Portal/port-sim-species-validation-20.RDS")
port.sim.sp30 <- readRDS("Scripts/HMMs/Portal/port-sim-species-validation-30.RDS")

port.sim.sp <- rbind(port.sim.sp15, port.sim.sp20, port.sim.sp30)

port.sim.sp <- as.data.frame(port.sim.sp)

colnames(port.sim.sp)[2] <- "old.code"

port.sim.sp <- merge(port.sim.sp, traits[, c(1,2)], by.x = "old.code", by.y = "species", all.x = T, all.y = F)

port.sim.sp <- port.sim.sp %>% 
  dplyr::mutate(across(sim:iter, ~as.numeric(.x)))

colnames(port.sim.sp)[4:9] <- paste(colnames(port.sim.sp)[4:9], "est", sep=".")

port.sim.sp <- merge(port.sim.sp, boot[,c(1:7)], by = "Code")

port.sim.sp <- merge(port.sim.sp, port.df[,c(2,10:13)], by = "Code")

port.rmse <- port.sim.sp %>%
  dplyr::group_by(sim, patches) %>%
  dplyr::filter(mean.years > 1) %>%
  dplyr::summarize(s.bias = mean(abs(s.boot - mean(s.est))),
            s.rmse = rmse(s.boot, s.est),
            s.var = var(s.boot, s.est),
            c.bias = mean(abs(c.boot - mean(c.est))),
            c.rmse = rmse(c.boot, c.est),
            c.var = var(c.boot, g.est),
            g.bias = mean(abs(g.boot - mean(g.est))),
            g.rmse = rmse(g.boot, g.est),
            g.var = var(g.boot, g.est)
            )

a <- ggplot(port.rmse, aes(x = s.rmse, group = patches, col = patches)) +
  geom_density(adjust = 1.5) +
  geom_vline(xintercept = mean(port.rmse[port.rmse$patches == 15,]$s.rmse), col = '#1b9e77', linetype = 2) +
  geom_vline(xintercept = mean(port.rmse[port.rmse$patches == 20,]$s.rmse), col = '#d95f02', linetype = 2) +
  geom_vline(xintercept = mean(port.rmse[port.rmse$patches == 30,]$s.rmse), col = '#7570b3', linetype = 2) +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) + theme_classic() +
  theme(
      legend.position = c(0.2,0.8)
  ) +
  labs(x = "RMSE of temporal dispersal estimates") +
  xlim(0,0.4)

b <- ggplot(port.rmse, aes(x = c.rmse, group = patches, col = patches)) +
  geom_density(adjust = 1.5) +
  geom_vline(xintercept = mean(port.rmse[port.rmse$patches == 15,]$c.rmse), col = '#1b9e77', linetype = 2) +
  geom_vline(xintercept = mean(port.rmse[port.rmse$patches == 20,]$c.rmse), col = '#d95f02', linetype = 2) +
  geom_vline(xintercept = mean(port.rmse[port.rmse$patches == 30,]$c.rmse), col = '#7570b3', linetype = 2) +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) + theme_classic() +
  theme(
      legend.position = c(0.8,0.8)
  ) +
  labs(x = "RMSE of spatial dispersal estimates") +
  xlim(0,0.3)

port.sims.fig <- ggarrange(a,b)

port.sims.fig <- annotate_figure(port.sims.fig, top = text_grob("Portal (max 64 patches, 22 years)", color = "black", face = "bold", size = 16))

ggplot(port.rmse, aes(x = s.bias, group = patches, col = patches)) +
  geom_density(adjust = 1.5)

ggplot(port.rmse, aes(x = s.var, group = patches, col = patches)) +
  geom_density(adjust = 1.5)

#### Carrizo Plain ####
##### Boots ####
cp.boot <- readRDS("Scripts/HMMs/Carrizo-Plain/CP-species-boot.RDS")
CP.meta <- read.csv("Data/Long-Term-Datasets/Carrizo-Plain/HMM-meta-CP.csv")

boot <- as.data.frame(cp.boot)

boot <- boot %>% 
  dplyr::mutate(across(p0.boot:perc.con, ~as.numeric(.x)))

colnames(boot)[1] <- "Code"

boot <- merge(boot, CP.df, by = "Code")

boot <- boot[boot$s.boot.sd < 0.1 & boot$c.boot.sd < 0.1 & boot$perc.con > 0.5,]

ggplot(boot[boot$mean.years > 1 & boot$n.plots >= 20,], aes(x = s.boot, y = c.boot)) +
  geom_smooth(method = "lm") +
  geom_text(aes(label = Code)) +
  geom_errorbar(aes(xmin = s.boot - s.boot.sd, xmax = s.boot + s.boot.sd)) +
  geom_errorbar(aes(ymin = c.boot - c.boot.sd, ymax = c.boot + c.boot.sd))

##### Patches #####
CP.sim.sp15 <- readRDS("Scripts/HMMs/Carrizo-Plain/CP-sim-species-validation-15.RDS")
CP.sim.sp20 <- readRDS("Scripts/HMMs/Carrizo-Plain/CP-sim-species-validation-20.RDS")
CP.sim.sp30 <- readRDS("Scripts/HMMs/Carrizo-Plain/CP-sim-species-validation-30.RDS")

CP.sim.sp <- rbind(CP.sim.sp15, CP.sim.sp20, CP.sim.sp30)

CP.sim.sp <- as.data.frame(CP.sim.sp)

colnames(CP.sim.sp)[2] <- "Code"

CP.sim.sp <- CP.sim.sp %>% 
  dplyr::mutate(across(sim:iter, ~as.numeric(.x)))

colnames(CP.sim.sp)[4:9] <- paste(colnames(CP.sim.sp)[4:9], "est", sep=".")

CP.sim.sp <- merge(CP.sim.sp, boot[,1:7], by = "Code")

CP.sim.sp <- merge(CP.sim.sp, CP.df[,c(2, 10:13)], by = "Code")

CP.rmse <- CP.sim.sp %>%
  dplyr::group_by(sim, patches) %>%
  dplyr::filter(mean.years > 1) %>%
  dplyr::summarize(s.bias = mean(abs(s.boot - mean(s.est))),
            s.rmse = rmse(s.boot, s.est),
            s.var = var(s.boot, s.est),
            c.bias = mean(abs(c.boot - mean(c.est))),
            c.rmse = rmse(c.boot, c.est),
            c.var = var(c.boot, g.est),
            g.bias = mean(abs(g.boot - mean(g.est))),
            g.rmse = rmse(g.boot, g.est),
            g.var = var(g.boot, g.est)
            )

a <- ggplot(CP.rmse, aes(x = s.rmse, group = patches, col = patches)) +
  geom_density(adjust = 1.5) +
  geom_vline(xintercept = mean(CP.rmse[CP.rmse$patches == 15,]$s.rmse), col = '#1b9e77', linetype = 2) +
  geom_vline(xintercept = mean(CP.rmse[CP.rmse$patches == 20,]$s.rmse), col = '#d95f02', linetype = 2) +
  geom_vline(xintercept = mean(CP.rmse[CP.rmse$patches == 30,]$s.rmse), col = '#7570b3', linetype = 2) +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) + theme_classic() +
  theme(
      legend.position = c(0.2,0.8)
  ) +
  labs(x = "RMSE of temporal dispersal estimates") +
  xlim(0,0.4)

b <- ggplot(CP.rmse, aes(x = c.rmse, group = patches, col = patches)) +
  geom_density(adjust = 1.5) +
  geom_vline(xintercept = mean(CP.rmse[CP.rmse$patches == 15,]$c.rmse), col = '#1b9e77', linetype = 2) +
  geom_vline(xintercept = mean(CP.rmse[CP.rmse$patches == 20,]$c.rmse), col = '#d95f02', linetype = 2) +
  geom_vline(xintercept = mean(CP.rmse[CP.rmse$patches == 30,]$c.rmse), col = '#7570b3', linetype = 2) +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3')) + theme_classic() +
  theme(
      legend.position = c(0.2,0.8)
  ) +
  labs(x = "RMSE of spatial dispersal estimates") +
  xlim(0,0.3)

CP.sims.fig <- ggarrange(a,b)

CP.sims.fig <- annotate_figure(CP.sims.fig, top = text_grob("Carrizo Plain (max 36 patches, 14 years)", color = "black", face = "bold", size = 16))

#### Jasper Ridge ####
##### Boots #####
boot <- readRDS("Scripts/HMMs/Jasper-Ridge/JR-species-boot-A.RDS")
boot2 <- readRDS("Scripts/HMMs/Jasper-Ridge/JR-species-boot-B.RDS")
trait <- read.csv("Data/Long-Term-Datasets/Jasper-Ridge/JR_Traits.csv")
meta.jrA <- read.csv("Data/Long-Term-Datasets/Jasper-Ridge/HMM-meta-JR-A.csv")
meta.jrB <- read.csv("Data/Long-Term-Datasets/Jasper-Ridge/HMM-meta-JR-B.csv")

meta.jrA$control <- "A"
meta.jrB$control <- "B"

meta.jr <- rbind(meta.jrA, meta.jrB)

meta.jr <- meta.jr %>%
  group_by(species) %>%
  summarize(across(n.plots:max.years, ~mean(.x)))

boot <- as.data.frame(boot)
boot2 <- as.data.frame(boot2)

boot$control <- "A"
boot2$control <- "B"

boot <- rbind(boot, boot2)

boot <- boot %>% 
  dplyr::mutate(across(p0.boot:perc.con, ~as.numeric(.x)))

trait$Code.4 <- toupper(trait$Code.4)
colnames(boot)[1] <- "Code.4"
colnames(meta.jr)[1] <- "Code.4"

boot <- merge(boot, trait[,c(2,3)], by = "Code.4")
meta.jr <- merge(meta.jr, trait[,c(2,3)], by = "Code.4")

boot  <- boot[,c(16,15,2:14)]

boot <- boot[boot$s.boot.sd < 0.1 & boot$c.boot.sd < 0.1 & boot$perc.con > 0.5,]

boot <- filter(boot, s.boot.sd < 0.1, c.boot.sd < 0.1, perc.con > 0.5)

boot <- merge(boot, meta.jr[,2:6], by = "Code")

ggplot(boot[boot$mean.years > 1 & boot$n.plots >= 15,], aes(x = s.boot, y = c.boot, col = control, group = control)) +
  geom_smooth(method = "lm") +
  geom_text(aes(label = Code)) +
  geom_errorbar(aes(xmin = s.boot - s.boot.sd, xmax = s.boot + s.boot.sd)) +
  geom_errorbar(aes(ymin = c.boot - c.boot.sd, ymax = c.boot + c.boot.sd))

boot <- boot %>%
  group_by(Code) %>%
  summarize(across(p0.boot:r.boot, ~mean(.x)))

#saveRDS(boot, "Scripts/HMMs/Jasper-Ridge/JR-species-boot-avg.RDS")

##### Patches #####
JR.sim.sp10 <- readRDS("Scripts/HMMs/Jasper-Ridge/JR-sim-species-validation-10.RDS")
JR.sim.sp15 <- readRDS("Scripts/HMMs/Jasper-Ridge/JR-sim-species-validation-15.RDS")
JR.sim.sp20 <- readRDS("Scripts/HMMs/Jasper-Ridge/JR-sim-species-validation-20.RDS")
JR.sim.sp30 <- readRDS("Scripts/HMMs/Jasper-Ridge/JR-sim-species-validation-30.RDS")

JR.sim.sp <- rbind(JR.sim.sp10, JR.sim.sp15, JR.sim.sp20, JR.sim.sp30)

JR.sim.sp <- as.data.frame(JR.sim.sp)


JR.sim.sp <- JR.sim.sp %>% 
  dplyr::mutate(across(sim:iter, ~as.numeric(.x)))

colnames(JR.sim.sp)[4:9] <- paste(colnames(JR.sim.sp)[4:9], "est", sep=".")

colnames(JR.sim.sp)[2] <- "Code"

JR.sim.sp <- merge(JR.sim.sp, boot[,1:6], by = "Code", all.x = T, all.y = F)

JR.sim.sp <- merge(JR.sim.sp, meta.jr[,2:6], by = "Code")

JR.rmse <- JR.sim.sp %>%
  dplyr::group_by(sim, patches) %>%
  dplyr::filter(mean.years > 1) %>%
  dplyr::summarize(s.bias = mean(abs(s.boot - mean(s.est))),
            s.rmse = rmse(s.boot, s.est),
            s.var = var(s.boot, s.est),
            c.bias = mean(abs(c.boot - mean(c.est))),
            c.rmse = rmse(c.boot, c.est),
            c.var = var(c.boot, g.est),
            g.bias = mean(abs(g.boot - mean(g.est))),
            g.rmse = rmse(g.boot, g.est),
            g.var = var(g.boot, g.est)
            )

      
a <- ggplot(JR.rmse, aes(x = s.rmse, group = patches, col = patches)) +
  geom_density(adjust = 1.5) +
  geom_vline(xintercept = mean(JR.rmse[JR.rmse$patches == 10,]$s.rmse), col = '#e7298a', linetype = 2) +
  geom_vline(xintercept = mean(JR.rmse[JR.rmse$patches == 15,]$s.rmse), col = '#1b9e77', linetype = 2) +
  geom_vline(xintercept = mean(JR.rmse[JR.rmse$patches == 20,]$s.rmse), col = '#d95f02', linetype = 2) +
  geom_vline(xintercept = mean(JR.rmse[JR.rmse$patches == 30,]$s.rmse), col = '#7570b3', linetype = 2) +
  scale_color_manual(values = c('#e7298a', '#1b9e77','#d95f02','#7570b3')) + theme_classic() +
  theme(
      legend.position = c(0.8,0.8)
  ) +
  labs(x = "RMSE of temporal dispersal estimates") +
  xlim(0,0.4)

b <- ggplot(JR.rmse, aes(x = c.rmse, group = patches, col = patches)) +
  geom_density(adjust = 1.5) +
  geom_vline(xintercept = mean(JR.rmse[JR.rmse$patches == 10,]$c.rmse), col = '#e7298a', linetype = 2) +
  geom_vline(xintercept = mean(JR.rmse[JR.rmse$patches == 15,]$c.rmse), col = '#1b9e77', linetype = 2) +
  geom_vline(xintercept = mean(JR.rmse[JR.rmse$patches == 20,]$c.rmse), col = '#d95f02', linetype = 2) +
  geom_vline(xintercept = mean(JR.rmse[JR.rmse$patches == 30,]$c.rmse), col = '#7570b3', linetype = 2) +
  scale_color_manual(values = c('#e7298a', '#1b9e77','#d95f02','#7570b3')) + theme_classic() +
  theme(
      legend.position = c(0.8,0.8)
  ) +
  labs(x = "RMSE of spatial dispersal estimates") +
  xlim(0,0.3)

JR.sims.fig <- ggarrange(a,b)

JR.sims.fig <- annotate_figure(JR.sims.fig, top = text_grob("Jasper Ridge (max 36 patches, 37 yrs)", color = "black", face = "bold", size = 16))

#### McLaughlin ####
##### Boots #####
boot <- readRDS("Scripts/HMMs/McLaughlin/McL-species-boot.RDS")
mcl.meta <- read.csv("Data/Long-Term-Datasets/McLaughlin/HMM-meta-mcl.csv")
boot <- as.data.frame(boot)

boot <- boot %>% 
  dplyr::mutate(across(p0.boot:perc.con, ~as.numeric(.x)))

boot <- boot[boot$s.boot.sd < 0.1 & boot$c.boot.sd < 0.1 & boot$perc.con > 0.5,]

boot <- merge(boot, mcl.df[,c(1,8:12)], by = "Species_Name")

ggplot(boot[boot$mean.years > 2 & boot$n.plots >= 20,], aes(x = s.boot, y = c.boot)) +
  geom_smooth(method = "lm") +
  geom_text(aes(label = Code)) +
  geom_errorbar(aes(xmin = s.boot - s.boot.sd, xmax = s.boot + s.boot.sd)) +
  geom_errorbar(aes(ymin = c.boot - c.boot.sd, ymax = c.boot + c.boot.sd))

##### Patches #####
mcl.sim.sp10 <- readRDS("Scripts/HMMs/McLaughlin/mcl-sim-species-validation-10.RDS")
mcl.sim.sp15 <- readRDS("Scripts/HMMs/McLaughlin/mcl-sim-species-validation-15.RDS")
mcl.sim.sp20 <- readRDS("Scripts/HMMs/McLaughlin/mcl-sim-species-validation-20.RDS")
mcl.sim.sp30 <- readRDS("Scripts/HMMs/McLaughlin/mcl-sim-species-validation-30.RDS")
mcl.sim.sp50 <- readRDS("Scripts/HMMs/McLaughlin/mcl-sim-species-validation-50.RDS")

mcl.sim.sp <- rbind(mcl.sim.sp10, mcl.sim.sp15, mcl.sim.sp20, mcl.sim.sp30, mcl.sim.sp50)

mcl.sim.sp <- as.data.frame(mcl.sim.sp)

mcl.sim.sp <- mcl.sim.sp %>% 
  dplyr::mutate(across(sim:iter, ~as.numeric(.x)))

colnames(mcl.sim.sp)[4:9] <- paste(colnames(mcl.sim.sp)[4:9], "est", sep=".")

mcl.sim.sp <- merge(mcl.sim.sp, boot[,c(1:7,15)], by = "Species_Name")

mcl.sim.sp <- merge(mcl.sim.sp, mcl.df[,c(1,8:12)], by = "Species_Name")

mcl.rmse <- mcl.sim.sp %>%
  dplyr::group_by(sim, patches) %>%
  dplyr::filter(mean.years > 1) %>%
  dplyr::summarize(s.bias = mean(abs(s.boot - mean(s.est))),
            s.rmse = rmse(s.boot, s.est),
            s.var = var(s.boot, s.est),
            c.bias = mean(abs(c.boot - mean(c.est))),
            c.rmse = rmse(c.boot, c.est),
            c.var = var(c.boot, g.est),
            g.bias = mean(abs(g.boot - mean(g.est))),
            g.rmse = rmse(g.boot, g.est),
            g.var = var(g.boot, g.est)
            )

a <- ggplot(mcl.rmse, aes(x = s.rmse, group = patches, col = patches)) +
  geom_density(adjust = 1.5) +
  geom_vline(xintercept = mean(mcl.rmse[mcl.rmse$patches == 10,]$s.rmse), col = '#e7298a', linetype = 2) +
  geom_vline(xintercept = mean(mcl.rmse[mcl.rmse$patches == 15,]$s.rmse), col = '#1b9e77', linetype = 2) +
  geom_vline(xintercept = mean(mcl.rmse[mcl.rmse$patches == 20,]$s.rmse), col = '#d95f02', linetype = 2) +
  geom_vline(xintercept = mean(mcl.rmse[mcl.rmse$patches == 30,]$s.rmse), col = '#7570b3', linetype = 2) +
  scale_color_manual(values = c('#e7298a', '#1b9e77','#d95f02','#7570b3')) + theme_classic() +
  theme(
      legend.position = c(0.8,0.8)
  ) +
  labs(x = "RMSE of temporal dispersal estimates") +
  xlim(0,0.4)

b <- ggplot(mcl.rmse, aes(x = c.rmse, group = patches, col = patches)) +
  geom_density(adjust = 1.5) +
  geom_vline(xintercept = mean(mcl.rmse[mcl.rmse$patches == 10,]$c.rmse), col = '#e7298a', linetype = 2) +
  geom_vline(xintercept = mean(mcl.rmse[mcl.rmse$patches == 15,]$c.rmse), col = '#1b9e77', linetype = 2) +
  geom_vline(xintercept = mean(mcl.rmse[mcl.rmse$patches == 20,]$c.rmse), col = '#d95f02', linetype = 2) +
  geom_vline(xintercept = mean(mcl.rmse[mcl.rmse$patches == 30,]$c.rmse), col = '#7570b3', linetype = 2) +
  scale_color_manual(values = c('#e7298a', '#1b9e77','#d95f02','#7570b3')) + theme_classic() +
  theme(
      legend.position = c(0.8,0.8)
  ) +
  labs(x = "RMSE of spatial dispersal estimates") +
  xlim(0,0.3)

mcl.sims.fig <- ggarrange(a,b)

mcl.sims.fig <- annotate_figure(mcl.sims.fig, top = text_grob("McLaughlin (max 400 patches, 19 years)", color = "black", face = "bold", size = 16))

ggplot(mcl.rmse, aes(x = s.bias, group = patches, col = patches)) +
  geom_density(adjust = 1.5)

ggplot(mcl.rmse, aes(x = s.var, group = patches, col = patches)) +
  geom_density(adjust = 1.5)

#### All Patches Figs ####
sim.fig <- ggarrange(SD.sims.fig, port.sims.fig, CP.sims.fig, JR.sims.fig, mcl.sims.fig, nrow = 5)

ggsave("Scripts/HMMs/Validation/sim-patches.jpg", sim.fig, height = 18, width = 8, units = "in", dpi = 600)
