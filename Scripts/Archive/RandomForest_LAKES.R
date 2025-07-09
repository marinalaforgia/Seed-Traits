#rm(list=ls())
########### LIBRARIES ###########
library(tidyverse)
library(rfUtilities)
library(randomForest)
library(caret)
library(vegan)
library(ggpubr)
#library(jtools)
library(corrplot)
cbPalette <- c("#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#E69F00", "#CC79A7", "#999999")

DATA <- Richness_complete

### REGULAR RANDOM FOREST
#Put your response variable name in here
Response<-"Rich"
# DATA$Dataset <- as.factor(DATA$Dataset)
#Put a dataframe with only index covariates in here 
Covariates<- DATA[,14:42]

tmp <- DATA %>% select(ACI, BKdB_low, BKdB_bird, H, SoundScapeI, ADI, AR, wind, avgAMP)
corrplot(cor(tmp), method = "number")

summary(lm(Rich~Roughness, data=DATA))
summary(lm(Rich~ACI, data=DATA))

RandomForestMod<-randomForest(as.formula(paste(Response, "~", paste(colnames(tmp), collapse = " + "),sep = "")), data=DATA, importance=TRUE, ntree=500)


### TAKE OUT CORRELATED INDEX COVARIATES
# using pearson correlations
tmp <- cor(Covariates)
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0
data.new <- Covariates[,!apply(tmp,2,function(x) any(abs(x) > 0.9))]


p_mat <- cor_5$p
corrplot(cor(data.new), type = "upper", order = "hclust", method="circle")

cl<-multi.collinear(DATA[,14:42], p=0.1) # have to exclude Dataset from the model bc its a factor
DATAcomplete <- DATA[,-which(colnames(DATA) %in% cl)] 

#Put a dataframe with only index covariates that ARE NOT collinear {still includes Dataset as a predictor}
Covariates2<- DATAcomplete[,14:34] 

# does not contain colinear variables, DOES contain DATASET
RandomForestMod2<-randomForest(as.formula(paste(Response, "~", paste(colnames(Covariates2), collapse = " + "),sep = "")), data=DATAcomplete, importance=TRUE, ntree=500)

#### RANDOM FOREST MODEL SELECTION
RFMODSEL<-rf.modelSel(DATAcomplete[,which(colnames(DATAcomplete)=="ACI"):which(colnames(DATAcomplete)=="AR")], DATAcomplete[[Response]],  imp.scale = "mir", final.model = TRUE)


##Get R squared and Mean squared error for each RF type
Output<-rbind(data.frame(MSE=mean((predict(RandomForestMod)-DATA[[Response]])^2), 
                         Rsq=mean(RandomForestMod$rsq), Model="RandomForest-Global", ResponseType=Response))
Output2<-rbind(data.frame(MSE=mean((predict(RandomForestMod2)-DATAcomplete[[Response]])^2), 
                         Rsq=mean(RandomForestMod2$rsq), Model="RandomForest-Global", ResponseType=Response))

Output_RFSelection<-rbind(data.frame(MSE=mean((RFMODSEL$rf.final$mse)), 
                                     Rsq=mean(RFMODSEL$rf.final$rsq), Model="RandomForest-Model Selection", ResponseType=Response, Dataset="Lakes"))



# this below is RB's code, I get an error about object of longer length is not multiple of shorter length
# Output_RFSelection<-data.frame(MSE=(RFMODSEL$rf.final$mse-DATAcomplete[[Response]]^2), 
#                    Rsq=RFMODSEL$rf.final$rsq, Model="RandomForest_ModelSelection", ResponseType=Response)

##Make predictions from each model type for your entire data set
#PredictionsSR_RFMODSEL<-predict(RFMODSEL$rf.final, newdata=Entiredataset)
#PredictionsSR_RFMOD<-predict(RandomForestMod, newdata=Entiredataset)

##Importance of each index in the predictive model

#Quick way
varimp.modsel <- varImp(RFMODSEL$rf.final)

varimp.modsel$Index <- rownames(varimp.modsel)
#varImp(RandomForestMod)
varImp(RandomForestMod2)

# for selected model
ggplot(varimp.modsel) +
  geom_point(aes(x=reorder(Index, Overall),y=Overall)) +
  coord_flip() +
  theme_minimal() +
  #geom_hline(aes(yintercept=10, color="red")) +
  labs(y="% increase in MSE if variable removed", x = "Acoustic Index", title="LAKES: Richness", subtitle = paste("Selected Model MSE=",Output_RFSelection$MSE,"; R^2=",Output_RFSelection$Rsq)) 
#title = "Relative importance of acoustic indices in predictive power of Final Random Forest")

######### observed vs predicted Species Richness plots #########


DATAcomplete$Rich_predict <- predict(RFMODSEL$rf.final)

DATAcomplete$Residuals <- DATAcomplete$Rich - DATAcomplete$Rich_predict
hist(DATAcomplete$Residuals) 

actuals_preds <- data.frame(cbind(actuals=DATAcomplete$Rich, predicteds=DATAcomplete$Rich_predict))  # make actuals_predicteds dataframe.
correlation_accuracy <- cor(actuals_preds) # 50% accurate ... model underpredicts richness
head(actuals_preds)

# relationship of Richness to each acoustic variable: many plots
varlist <- RFMODSEL$selvars

vars <- DATAcomplete %>% select(Rich, all_of(varlist))

for (i in varlist){
  plot <- ggplot(vars, aes_string(vars$Rich, x=i))+ 
    geom_point() +
    ylab("Richness") +
    theme_classic()
  print(plot)
}


ggplot(DATAcomplete) +
  geom_boxplot(stat="boxplot", aes(group=Rich, x=Rich,y=Residuals))

ggplot(DATAcomplete) +
  geom_point(aes(x=BKdB_bird,y=Residuals))

ggplot(DATAcomplete) +
  geom_point(aes(x=BKdB_bird,y=Rich))

ggplot(DATAcomplete) +
  geom_boxplot(stat="boxplot", aes(group=Site, x=Site,y=avgAMP)) 

ggplot(DATAcomplete,aes(Rich, Rich_predict)) +
  geom_point(,size=3, alpha=0.5) +
  labs(x="Observed Species Richness", y="Predicted Species Richness") +
  theme(aspect.ratio=1, legend.position="top", legend.text = element_text(size = 10), legend.background = element_blank()) +
  theme_classic() +
  scale_color_manual(values=cbPalette[1:3]) +
  scale_y_continuous(limits = c(0,17), breaks=seq(0,17,3), expand=c(0,0)) +
  scale_x_continuous(limits = c(0,17), breaks=seq(0,17,3), expand=c(0,0)) +
  geom_abline() 
#stat_smooth(aes(x=Richness, y=Rich_predict), formula = y ~ log2(x), method="lm")


######## Cross Validate each model {RF code} ##########

COVARIATES <- Covariates2

CrossValidationFunction<-function(DATAcomplete, Response, COVARIATES, p=.10, n=99){
  
  ##Get random subsamples of the dataframe
  Datawitheld<-round(p*nrow(DATAcomplete),0)
  Rownums<- replicate(n, sample(nrow(DATAcomplete), Datawitheld), simplify=FALSE)
  
  ##Loop through each subset as training data
  MSElist<-lapply(1:n, function(i){
    
    #Use test and train data in models
    testData<-DATAcomplete[Rownums[[i]], ]
    trainData<-DATAcomplete[-Rownums[[i]],]
    Test2<-randomForest(as.formula(paste(Response, "~", 
                                         paste(COVARIATES, collapse = " + "),sep = "")), data=trainData, 
                        importance=TRUE)
    
    Predictions<-predict(Test2, testData)
    MSE=mean((Predictions-testData[[Response]])^2)
    data.frame(MeanSqError=MSE, Iteration=i)
  })##End of loop through a priori models
  ##Bind altogether
  MSElist<-do.call("rbind", MSElist)
  MSElist
}


RFMODSELCV<-CrossValidationFunction(DATAcomplete, Response, RFMODSEL$selvars, p=.10, n=99)

summary(RFMODSELCV)
