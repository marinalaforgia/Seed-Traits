## Transforming oil content of KEW ##

sid <- read.csv("Protocols/Seed-traits/Datasets/Kew-SID_Oil_Full.csv")
sid$full.species <- paste(sid$Genus, sid$Species, sid$var.ssp, sid$var.ssp.name)
sid$full.species <- trimws(sid$full.species, which = "right")

sid.cbg<- sid[sid$full.species %in% cbg$name,]
sid.ca <- sid[sid$full.species %in% Calflora.DA$Species,]
sid.az <- sid[sid$full.species %in% AZ.annuals$Scientific.Name,]
sid.az <- sid.az[!sid.az$full.species %in% sid.ca$full.species,]

#write.csv(sid.ca, "Protocols/Seed-traits/Datasets/sid-ca.csv", row.names = F)
#write.csv(sid.az, "Protocols/Seed-traits/Datasets/sid-az.csv", row.names = F)

try <- read.csv("Protocols/Seed-traits/Datasets/TRY-Oil.csv")
try$Species <- paste(try$Genus, try$Species, try$var_ssp, try$vs_name)
try$Species <- trimws(try$Species, which = "right")

try.ca <- try[try$Species %in% Calflora.DA$Species,]
try.az <- try[try$Species %in% AZ.annuals$Scientific.Name,]
try.cbg<- try[try$Species %in% cbg$name,]

try.ps<- try[try$Species %in% pos.spp$Species,]
sid.ps<- sid[sid$full.species %in% pos.spp$Species,]

sid <- read.csv("Protocols/Seed-traits/Datasets/KEW-SID_Oil-Content_orig.csv")
write.csv(sid, "Protocols/Seed-traits/Datasets/sid.csv", row.names = F)

ipc <- read.csv("Protocols/Seed-traits/Datasets/Cal-IPC_Annual.csv")
sid.ipc <- sid[sid$full.species %in% ipc$Species,]
