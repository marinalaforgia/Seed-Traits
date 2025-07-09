#### Seed settling time cleaning ####
rm(list=ls())
set <- read.csv("Data/Settling-time/Settling-time_entered/20220805_Settling-time_entered.csv")

acc <- read.csv("Data/20221216_Seeds_All-Accessions.csv")
species <- read.csv("Data/20211001_Full-Species-List.csv")

set[set$Species == "Antirrhinum coulterianum" & set$Replicate == 2,]$Time <- 1.620 # this is an entry error that has been double checked

set$set.time.mpsec <- 185.42/set$Time # Tube is 73 inches = 185.42 cm, Kate double checked
set$Species <- trimws(set$Species, which = "both")


# brief cleaning
unique(set[!set$Species %in% species$Species,]$Species)

set$Species <- recode(set$Species, "Centuarea solstitialis" = "Centaurea solstitialis",
                              "Chenopodium berlandieri var. zschackii" = "Chenopodium berlandieri var. zschackei",
                              "Chylismia claviformis spp claviformis" = "Chylismia claviformis ssp. claviformis",
                              "Chylismia claviformis ssp funerea" = "Chylismia claviformis ssp. funerea",      
                              "Helianthus niveus ssp canescens" = "Helianthus niveus ssp. canescens",        
                              "Helianthus petiolaris ssp fallax" = "Helianthus petiolaris ssp. fallax",        
                              "Heliomeris longifolia var annua" =  "Heliomeris longifolia var. annua",         
                              "Layia glandolara" = "Layia glandulosa",                         
                              "Dysphania ambrosiodes" = "Dysphania ambrosioides",                   
                              "Eremothera boothi" = "Eremothera boothii",                        
                              "Eriastum eremicum" = "Eriastrum eremicum",                       
                              "Eriastum sapphirium" = "Eriastrum sapphirinum",                     
                              "Eriodium botrys" = "Erodium botrys",                          
                              "Claytonia parviflora ssp parviflora" = "Claytonia parviflora ssp. parviflora",      
                              "Claytonia parviflora ssp viridis" =  "Claytonia parviflora ssp. viridis",       
                              "Coroylanthus hevinii" = "Cordylanthus nevinii",                     
                              "Coroylanthus rigidus ssp setiger" = "Cordylanthus rigidus ssp. setiger",        
                              "Cryptantha nevadensis var rigida" =  "Cryptantha nevadensis var. rigida",       
                              "Crypthantha pterocarya" =  "Cryptantha pterocarya",                  
                              "Daucus Pusillus" = "Daucus pusillus",                        
                              "Lepedium nitidum" = "Lepidium nitidum",                        
                              "Microseris gracilis" = "Microsteris gracilis",                      
                              "Monarda citriodora var austromontana" = "Monarda citriodora var. austromontana",     
                              "Monardella breweri ssp lanceolata" = "Monardella breweri ssp. lanceolata",       
                              "Muhlenbergia microsperms" = "Muhlenbergia microsperma",                 
                              "Nemacladus longiflorus var breviflorus" =  "Nemacladus longiflorus var. breviflorus",  
                              "Pectocarya pencillata" = "Pectocarya penicillata",                   
                              "Spergulana platensis" =  "Spergularia platensis",                    
                              "Trichostema austromontanum ssp. austromo" =  "Trichostema austromontanum ssp. austromontanum",
                              "Triofolium hirtum" = "Trifolium hirtum",                      
                              "Triofolium microcephalum" = "Trifolium microcephalum",                 
                              "Triofolium willdenovii" =  "Trifolium willdenovii")
               
unique(set[!set$Species %in% species$Species,]$Species)
unique(set[!set$Species %in% acc$Species,]$Species)

cultivated <- unique(acc[acc$site == "GRIN" & acc$type == "cultivated",]$Species)

set <- filter(set, !Species %in% cultivated, Species != "Poa secunda")

length(unique(set$Species))

write.csv(set, "Data/Settling-time/Settling-time_clean/20221216_Settling-Time_clean.csv", row.names = F)
