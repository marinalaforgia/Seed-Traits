# try data
library(plyr)
library(tidyverse)
#try.sp <- read.csv("Projects/Seed-Traits/Try/Try-Species-list.csv")

#species <- read.csv("Projects/Seed-Traits/Traits/Embryo-Morphology/Species-List.csv")
species <- read.csv("/Users/Marina/Documents/USDA-PostDoc/Projects/RGR-WUE/McL_80SitesSpeciesTraits_012615.csv")
phen <- read.csv("/Users/Marina/Documents/UC-Davis/Projects/Explorations/McL-Phenology/Data/McL_Phen_combined-barb-cal.csv")
species <- merge(species, phen[,1:2], by.x = "Species_Name_J12", by.y = "Species_Name")

#try.sp <- filter(try.sp, AccSpeciesName %in% species$Species)
#write.csv(try.sp, "Try-species-collection.csv", row.names = F)

try <- read.csv("/Users/Marina/Documents/USDA-PostDoc/Projects/Seed-Traits/Try/Try-Trait-Data.csv")

try <- filter(try, TraitID == 77 |  TraitID == 89 | TraitID == 3478) 
#try.2 <- merge(try, species, by.x = "SpeciesName", by.y = "Species_Name", all.x = F)
try <- merge(try, species, by.x = "SpeciesName", by.y = "Species_Name_J12", all.x = F)

try$StdValue <- ifelse(try$OriglName == "d13C", as.numeric(-0.008 - (try$StdValue * 0.001))/(1 + (try$StdValue * 0.001))*1000, as.numeric(try$StdValue))
try <- ddply(try, .(SpeciesName, TraitID), summarize, value = mean(StdValue))
try$trait <- ifelse(try$TraitID == 77, "RGR", "D13C")
try.D13C <- filter(try, trait == "D13C" )
try.D13C <- ddply(try.D13C, .(SpeciesName), summarize, DC13 = mean(value))



sla.try <- read.csv("/Users/Marina/Documents/USDA-PostDoc/Projects/RGR-WUE/sla-Dc13.csv")
sla.try <- filter(sla.try, Site != "TRY")

carrizo <- read.csv("/Users/Marina/Documents/USDA-PostDoc/Projects/RGR-WUE/Carrizo_plant trait data.csv")
carrizo <- filter(carrizo, ann.peren == "annual")
ggplot(carrizo, aes(y = log(SLA.cm2.g), x = D13C, col = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  scale_x_reverse() +
  geom_label(aes(label = SpeciesCode))

carrizo.pca <- carrizo[complete.cases(carrizo),]
pc.car <- prcomp(carrizo.pca[,c(2,4,5,6,8)], scale = T)
summary(pc.car)
biplot(pc.car)


autoplot(pc.car, data = carrizo.pca, loadings = T, loadings.label = T, label = T, shape = F, col = 'Status') +
  theme_classic() +
  scale_color_manual(values = c("magenta4","#1F968BFF", "goldenrod3"))

sla.all <- merge(carrizo)
sla.try$site.nat.inv <- paste(sla.try$Site, sla.try$nat.inv, sep = ".")

ggplot(sla.try[sla.try$System != "Carrizo",], aes(y = log(SLA), x = D13C, col = nat.inv, group = nat.inv)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  scale_x_reverse() +
  # geom_label(aes(label = Species.Short), size = 2) +
  # facet_wrap(~nat.inv) +
  labs(y = "Log of Specific Leaf Area", x = "Carbon Isotope Discrimination")

sla.try$D13C <- sla.try$D13C*-1

m1 <- lm(log(SLA) ~ D13C * Site, sla.try)
plot(fitted(m1), resid(m1))
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)

try.w <- pivot_wider(try, names_from = "trait", values_from = "value")
species <- merge(species, try.w[,c(1,4)], by.y = "SpeciesName", by.x = "Species_Name_J12", all.x = T)

try.MCL <- read.csv("/Users/Marina/Documents/USDA-PostDoc/Projects/RGR-WUE/try_McL.csv")

try.MCL <- merge(try.MCL, sla)
ggplot(try.MCL, aes(y = log(SLANS), x = D13C, group = Native.Exotic, col = Native.Exotic)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()


