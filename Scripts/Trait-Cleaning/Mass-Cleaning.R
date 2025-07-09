#### Mass Cleaning ####
rm(list=ls())
# For dimporphic seeds I will be using the average ###
# As of 08/04/2022: species still need to be done: THYCUR, MEDPOL, ANOPEN, CHOBREB, VOLTUB (chemical only)
# not going to use cultivated seeds from GRIN, or Poa (perennial)
# all species have a morphological unit, but only some species have a chemical unit; chemical unit just describes what level nutrient analysis are done on if this is DIFFERENT from the morphological unit (i.e., AVEBAR has a caryopsis inside a modified floral case, the full seed plus persistent case is the morphological unit, the caryopsis is the chemical unit)

acc <- read.csv("Data/20230323_Seeds_All-Accessions.csv", na.strings = "")
acc <- filter(acc, !is.na(use.))
mass <- read.csv("Data/Seed-mass/Seed-mass_entered/20230323_seed-mass_entered.csv", na.strings = "")

cultivated <- unique(acc[acc$site == "GRIN" & acc$type == "cultivated",]$Species)

species1 <- acc$Species[1:25]
species2 <- acc$Species[26:50]
species3 <- acc$Species[51:100]
species4 <- acc$Species[101:150]
species5 <- acc$Species[151:200]
species6 <- acc$Species[201:250]
species7 <- acc$Species[251:300]
species8 <- acc$Species[301:350]
species9 <- acc$Species[351:365]

ggplot(mass[mass$chem.morph.updated == "morph" & mass$Species %in% species9,], aes(y = Mass.mg, x = Species)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mass <- merge(mass, acc[,c(10,26,27)], by = "ID")

mass <- filter(mass, N > 5, !Species %in% cultivated, Species != "Poa secunda")

write.csv(mass, "Data/Seed-mass/Seed-mass_cleaned/20230323_seed-mass_cleaned.csv", row.names = F)
# all species have a morphological unit, but only some species have a chemical unit; chemical unit just describes what level nutrient analysis are done on if this is DIFFERENT from the morphological unit (i.e., AVEBAR has a caryopsis inside a modified floral case, the full seed plus persistent case is the morphological unit, the caryopsis is the chemical unit)


