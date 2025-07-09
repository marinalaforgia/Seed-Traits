traits <- read.csv("Data/20220805_Seed-Traits_cleaning.csv")

# PCA needs normalized data, these worked for my larger datasets, may be unnecessary or may work for just the CP data
traits$set.time.mpsec.sqrt <- sqrt(traits$set.time.mpsec)
traits$shape.asinsqrt <- asin(sqrt(traits$shape))
traits$mass.mg.log <- log(traits$mass.mg)
traits$height.cm.log <- log(traits$height.cm)
traits$wing.loading.log <- log(traits$wing.loading)
traits$cn.sqrt <- sqrt(traits$cn)
traits$C.prop.asinsqrt <- asin(sqrt(traits$prop.C))
traits$N.prop.log <- log(traits$prop.N)
traits$size.log <- log10(traits$size)

### PCA ####
trait.pca <- traits[complete.cases(traits[, -17]), c(1:5, 11,13, 17, 19:27)] # PCA needs complete data (no NAs), these columns will likely have to change
pca <- prcomp(trait.pca[, c(9:14,17)], scale = T) # choose columns to run PCA on
summary(pca) # maount of variance explained by each PC axis
biplot(pca) # biplot of pca

autoplot(pca, x = 1, y = 2, data = trait.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "appendage.type", shape = "appendage.type") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) # better visualization of data

trait.pca <- cbind(trait.pca, pca$x[,1:3]) # bind PC axes with trait data

# from here you could bind with your seed bank abundace dataset and look at correlations between the PCs and abundance or individual traits and abundance to see if there are any patterns