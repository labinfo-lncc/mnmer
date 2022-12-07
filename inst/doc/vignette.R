## ---- eval=FALSE--------------------------------------------------------------
#  library(devtools)
#  install_github("labinfo-lncc/mnmer", ref="main")

## ---- eval=TRUE---------------------------------------------------------------
library("mnmer")
dir <-system.file("extdata", package="mnmer")

## ---- eval=TRUE---------------------------------------------------------------
mosquito <- mnmer(file.path(dir, "mosquito_vir.fasta"),2,0)
plant <- mnmer(file.path(dir, "plant_vir.fasta"),2,0)

## ---- eval=TRUE---------------------------------------------------------------
mosquito <- mnmer(file.path(dir, "mosquito_vir.fasta"),2,1)
plant <- mnmer(file.path(dir, "plant_vir.fasta"),2,1)

## ---- eval=TRUE---------------------------------------------------------------
library(caret)
classes <- replicate(nrow(mosquito), "mosquito.vir")
featureMatrix_mosquito <- cbind(mosquito,classes)
classes <- replicate(nrow(plant), "plant.vir")
featureMatrix_plant <- cbind(plant,classes)

featureMatrix <- rbind(featureMatrix_mosquito, featureMatrix_plant)
featureMatrix <- subset(featureMatrix, select = -c(seqid))
train_index <- caret::createDataPartition(featureMatrix$classes, p=0.8, list=FALSE)
train <- featureMatrix[train_index, ]
test <- featureMatrix[-train_index, ]
control <- caret::trainControl(method="cv", 
                        summaryFunction=twoClassSummary, 
                        classProbs=TRUE, 
                        savePredictions = TRUE)
roc <- caret::train(classes ~ ., 
            data=train, 
            method="rf", 
            preProc=c("center"), 
            trControl=control)
res <- MLeval::evalm(roc) # Make the ROC plot

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

