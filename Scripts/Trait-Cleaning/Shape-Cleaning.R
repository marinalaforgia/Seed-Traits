#### Shape Cleaning ####
# Shape describes how round something is, the closer it is to 0 the more round it is, the closer it is to 0.3 the more needle like it is

rm(list=ls())

acc <- read.csv("Data/20230323_Seeds_All-Accessions.csv")
shape <- read.csv("Data/Shape_Size/20230322_Videometer.csv")

cultivated <- unique(acc[acc$site == "GRIN" & acc$type == "cultivated",]$Species)
shape <- filter(shape, !Species %in% cultivated, Species != "Poa secunda", Species != "Nicotiana attenuata")

shape.sum <- ddply(shape, .(Species, Species.Code, ID, chem.morph.fun, X2d.3d, app), summarize, area.mm2 = mean(Area..mm2.), length.mm = mean(Length..mm.), width.mm = mean(Width..mm.))

shape.sum.w <- pivot_wider(shape.sum, names_from = "X2d.3d", values_from = c("area.mm2", "length.mm", "width.mm"))

# Vs = Σ (xi – mean (x))2/n with n = 3 and x1 = length/length, x2 = height/length and x3 = width/length.
shape.sum.w$avg.length.mm <- ifelse(is.na(shape.sum.w$length.mm_3D), shape.sum.w$length.mm_2D, (shape.sum.w$length.mm_2D + shape.sum.w$length.mm_3D)/2)
shape.sum.w$thick.mm <- ifelse(is.na(shape.sum.w$width.mm_3D), shape.sum.w$width.mm_2D, shape.sum.w$width.mm_3D)

# length = longest axis = avg.length.mm, width = width.mm_2d, thickness/height = thick.mm
shape.sum.w$size.mm <- shape.sum.w$avg.length.mm
shape.sum.w$x1 <- shape.sum.w$avg.length.mm/shape.sum.w$avg.length.mm
shape.sum.w$x2 <- shape.sum.w$thick.mm/shape.sum.w$avg.length.mm
shape.sum.w$x3 <- shape.sum.w$width.mm_2D/shape.sum.w$avg.length.mm
shape.sum.w$x <- (shape.sum.w$x1 + shape.sum.w$x2 + shape.sum.w$x3)/3
shape.sum.w$Vs1 <- (shape.sum.w$x1 - shape.sum.w$x)^2 
shape.sum.w$Vs2 <- (shape.sum.w$x2 - shape.sum.w$x)^2 
shape.sum.w$Vs3 <- (shape.sum.w$x3 - shape.sum.w$x)^2 
shape.sum.w$shape <- (shape.sum.w$Vs1 + shape.sum.w$Vs2 + shape.sum.w$Vs3)/3
hist(shape.sum.w$shape)

shape.sum.w$full.area.mm2 <- ifelse(is.na(shape.sum.w$area.mm2_3D), shape.sum.w$area.mm2_2D*4, shape.sum.w$area.mm2_2D*2 + shape.sum.w$area.mm2_3D*2)

shape.sum.w <- shape.sum.w[,c(1:5, 13, 14, 22, 23)]

length(unique(shape.sum.w$Species))
colnames(shape.sum.w)[6] <- "width.mm"


family <- read.csv("Data/20211001_Full-Species-List.csv")

shape.sum.w <- merge(shape, family[,c(2,8)], by = "Species", all.x = T, all.y = F)

write.csv(shape.sum.w, "Data/Shape_Size/20230323_Shape-Size_Cleaned.csv", row.names = F)
