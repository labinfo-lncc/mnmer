library(data.table)
library(doParallel)
library(dplyr)
library(randomForest)
library(caret)
library(ranger)
library(e1071)

cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

et_grid =  expand.grid(mtry = 4:7, numRandomCuts = 1:10)
control <- trainControl(method="repeatedcv", number=10, repeats=3)
ctrl <- trainControl(method="cv", summaryFunction=twoClassSummary, classProbs=T, savePredictions = T)

#### Run classification from the best features (selected by RFE)
## Random Forest (ExtraTrees)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/Kmer3_21_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_genome_class2 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_genome_class2 <- train(classes ~ ., data=train, method="extraTrees", preProc=c("center", "scale"), trControl=ctrl,tuneGrid = et_grid)
saveRDS(model_genome_class2, "model_genome_class2_roc.rds")
pred <- predict(model_genome_class2, newdata = test)
cf <- confusionMatrix(pred, as.factor(test$classes))
byClass <- cf$byClass
overall <- cf$overall
VP <- cf$table[1,1]
FN <- cf$table[1,2]
FP <- cf$table[2,1]
VN <- cf$table[2,2]
table_tmp <- rbind(VP,FN)
table_tmp <- rbind(table_tmp,FP)
table_tmp <- rbind(table_tmp,VN)
table <- cbind(table, table_tmp)
write.table(overall, file="./metrics/Accuracy_Model_Genome_Class2.txt", row.names=F)
write.table(byClass, file="./metrics/byClass_Model_Genome_Class2.txt", row.names=F)
write.table(pred, file="./metrics/predictions_Model_Genome_Class2.txt", row.names=F)
write.table(table, file="./metrics/confusionMatrix_Model_Genome_Class2.txt", row.names=F)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/nonKmer4_13_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_genome_class1 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_genome_class1 <- train(classes ~ ., data=train, method="extraTrees", preProc=c("center", "scale"), trControl=ctrl,tuneGrid = et_grid)
saveRDS(model_genome_class1, "model_genome_class1_roc.rds")
pred <- predict(model_genome_class1, newdata = test)
cf <- confusionMatrix(pred, as.factor(test$classes))
byClass <- cf$byClass
overall <- cf$overall
VP <- cf$table[1,1]
FN <- cf$table[1,2]
FP <- cf$table[2,1]
VN <- cf$table[2,2]
table_tmp <- rbind(VP,FN)
table_tmp <- rbind(table_tmp,FP)
table_tmp <- rbind(table_tmp,VN)
table <- cbind(table, table_tmp)
write.table(overall, file="./metrics/Accuracy_Model_Genome_Class1.txt", row.names=F)
write.table(byClass, file="./metrics/byClass_Model_Genome_Class1.txt", row.names=F)
write.table(pred, file="./metrics/predictions_Model_Genome_Class1.txt", row.names=F)
write.table(table, file="./metrics/confusionMatrix_Model_Genome_Class1.txt", row.names=F)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/Kmer3_12_5000bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_5000bp_class2 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_5000bp_class2 <- train(classes ~ ., data=train, method="extraTrees", preProc=c("center", "scale"), trControl=ctrl,tuneGrid = et_grid)
saveRDS(model_5000bp_class2, "model_5000bp_class2_roc.rds")
pred <- predict(model_5000bp_class2, newdata = test) 
cf <- confusionMatrix(pred, as.factor(test$classes))
byClass <- cf$byClass
overall <- cf$overall
VP <- cf$table[1,1]
FN <- cf$table[1,2]
FP <- cf$table[2,1]
VN <- cf$table[2,2]
table_tmp <- rbind(VP,FN)
table_tmp <- rbind(table_tmp,FP)
table_tmp <- rbind(table_tmp,VN)
table <- cbind(table, table_tmp)
write.table(overall, file="./metrics/Accuracy_Model_5000bp_Class2.txt", row.names=F)
write.table(byClass, file="./metrics/byClass_Model_5000bp_Class2.txt", row.names=F)
write.table(pred, file="./metrics/predictions_Model_5000bp_Class2.txt", row.names=F)
write.table(table, file="./metrics/confusionMatrix_Model_5000bp_Class2.txt", row.names=F)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/nonKmer3_12_5000bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
#subseq <- subseq[,-1]
subseq <- subset(subseq, select = -c(V1))
col_5000bp_class1 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_5000bp_class1 <- train(classes ~ ., data=train, method="extraTrees", preProc=c("center", "scale"), trControl=ctrl,tuneGrid = et_grid)
saveRDS(model_5000bp_class1, "model_5000bp_class1_roc.rds")
pred <- predict(model_5000bp_class1, newdata = test)
cf <- confusionMatrix(pred, as.factor(test$classes))
byClass <- cf$byClass
overall <- cf$overall
VP <- cf$table[1,1]
FN <- cf$table[1,2]
FP <- cf$table[2,1]
VN <- cf$table[2,2]
table_tmp <- rbind(VP,FN)
table_tmp <- rbind(table_tmp,FP)
table_tmp <- rbind(table_tmp,VN)
table <- cbind(table, table_tmp)
write.table(overall, file="./metrics/Accuracy_Model_5000bp_Class1.txt", row.names=F)
write.table(byClass, file="./metrics/byClass_Model_5000bp_Class1.txt", row.names=F)
write.table(pred, file="./metrics/predictions_Model_5000bp_Class1.txt", row.names=F)
write.table(table, file="./metrics/confusionMatrix_Model_5000bp_Class1.txt", row.names=F)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/Kmer3_12_3000bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_3000bp_class2 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_3000bp_class2 <- train(classes ~ ., data=train, method="extraTrees", preProc=c("center", "scale"), trControl=ctrl,tuneGrid = et_grid)
saveRDS(model_3000bp_class2, "model_3000bp_class2_roc.rds")
pred <- predict(model_3000bp_class2, newdata = test)
cf <- confusionMatrix(pred, as.factor(test$classes))
byClass <- cf$byClass
overall <- cf$overall
VP <- cf$table[1,1]
FN <- cf$table[1,2]
FP <- cf$table[2,1]
VN <- cf$table[2,2]
table_tmp <- rbind(VP,FN)
table_tmp <- rbind(table_tmp,FP)
table_tmp <- rbind(table_tmp,VN)
table <- cbind(table, table_tmp)
write.table(overall, file="./metrics/Accuracy_Model_3000bp_Class2.txt", row.names=F)
write.table(byClass, file="./metrics/byClass_Model_3000bp_Class2.txt", row.names=F)
write.table(pred, file="./metrics/predictions_Model_3000bp_Class2.txt", row.names=F)
write.table(table, file="./metrics/confusionMatrix_Model_3000bp_Class2.txt", row.names=F)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/Kmer3_12_1000bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_1000bp_class2 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_1000bp_class2 <- train(classes ~ ., data=train, method="extraTrees", preProc=c("center", "scale"), trControl=ctrl,tuneGrid = et_grid)
saveRDS(model_1000bp_class2, "model_1000bp_class2_roc.rds")
pred <- predict(model_1000bp_class2, newdata = test)
cf <- confusionMatrix(pred, as.factor(test$classes))
byClass <- cf$byClass
overall <- cf$overall
VP <- cf$table[1,1]
FN <- cf$table[1,2]
FP <- cf$table[2,1]
VN <- cf$table[2,2]
table_tmp <- rbind(VP,FN)
table_tmp <- rbind(table_tmp,FP)
table_tmp <- rbind(table_tmp,VN)
table <- cbind(table, table_tmp)
write.table(overall, file="./metrics/Accuracy_Model_1000bp_Class2.txt", row.names=F)
write.table(byClass, file="./metrics/byClass_Model_1000bp_Class2.txt", row.names=F)
write.table(pred, file="./metrics/predictions_Model_1000bp_Class2.txt", row.names=F)
write.table(table, file="./metrics/confusionMatrix_Model_1000bp_Class2.txt", row.names=F)

## KNN

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/nonKmer3_12_3000bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_3000bp_class1 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_3000bp_class1 <- train(classes ~ ., data=train, method = "knn", trControl = control, preProcess = c("center","scale"), tuneLength = 20)
saveRDS(model_3000bp_class1, "model_3000bp_class1_roc.rds")
pred <- predict(model_3000bp_class1, newdata = test) 
cf <- confusionMatrix(pred, as.factor(test$classes))
byClass <- cf$byClass
overall <- cf$overall
VP <- cf$table[1,1]
FN <- cf$table[1,2]
FP <- cf$table[2,1]
VN <- cf$table[2,2]
table_tmp <- rbind(VP,FN)
table_tmp <- rbind(table_tmp,FP)
table_tmp <- rbind(table_tmp,VN)
table <- cbind(table, table_tmp)
write.table(overall, file="./metrics/Accuracy_Model_3000bp_Class1.txt", row.names=F)
write.table(byClass, file="./metrics/byClass_Model_3000bp_Class1.txt", row.names=F)
write.table(pred, file="./metrics/predictions_Model_3000bp_Class1.txt", row.names=F)
write.table(table, file="./metrics/confusionMatrix_Model_3000bp_Class1.txt", row.names=F)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/nonKmer4_13_1000bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_1000bp_class1 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_1000bp_class1 <- train(classes ~ ., data=train, method = "knn", trControl = control, preProcess = c("center","scale"), tuneLength = 20)
saveRDS(model_1000bp_class1, "model_1000bp_class1_roc.rds")
pred <- predict(model_1000bp_class1, newdata = test)
cf <- confusionMatrix(pred, as.factor(test$classes))
byClass <- cf$byClass
overall <- cf$overall
VP <- cf$table[1,1]
FN <- cf$table[1,2]
FP <- cf$table[2,1]
VN <- cf$table[2,2]
table_tmp <- rbind(VP,FN)
table_tmp <- rbind(table_tmp,FP)
table_tmp <- rbind(table_tmp,VN)
table <- cbind(table, table_tmp)
write.table(overall, file="./metrics/Accuracy_Model_1000bp_Class1.txt", row.names=F)
write.table(byClass, file="./metrics/byClass_Model_1000bp_Class1.txt", row.names=F)
write.table(pred, file="./metrics/predictions_Model_1000bp_Class1.txt", row.names=F)
write.table(table, file="./metrics/confusionMatrix_Model_1000bp_Class1.txt", row.names=F)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/nonKmer3_12_500bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_500bp_class1 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_500bp_class1 <- train(classes ~ ., data=train, method = "knn", trControl = control, preProcess = c("center","scale"), tuneLength = 20)
saveRDS(model_500bp_class1, "model_500bp_class1_roc.rds")
pred <- predict(model_500bp_class1, newdata = test) 
cf <- confusionMatrix(pred, as.factor(test$classes))
byClass <- cf$byClass
overall <- cf$overall
VP <- cf$table[1,1]
FN <- cf$table[1,2]
FP <- cf$table[2,1]
VN <- cf$table[2,2]
table_tmp <- rbind(VP,FN)
table_tmp <- rbind(table_tmp,FP)
table_tmp <- rbind(table_tmp,VN)
table <- cbind(table, table_tmp)
write.table(overall, file="./metrics/Accuracy_Model_500bp_Class1.txt", row.names=F)
write.table(byClass, file="./metrics/byClass_Model_500bp_Class1.txt", row.names=F)
write.table(pred, file="./metrics/predictions_Model_500bp_Class1.txt", row.names=F)
write.table(table, file="./metrics/confusionMatrix_Model_500bp_Class1.txt", row.names=F)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/Kmer3_12_500bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
col_500bp_class2 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_500bp_class2 <- train(classes ~ ., data=train, method = "knn", trControl = control, preProcess = c("center","scale"), tuneLength = 20)
saveRDS(model_500bp_class2, "model_500bp_class2_roc.rds")
pred <- predict(model_500bp_class2, newdata = test)
cf <- confusionMatrix(pred, as.factor(test$classes))
byClass <- cf$byClass
overall <- cf$overall
VP <- cf$table[1,1]
FN <- cf$table[1,2]
FP <- cf$table[2,1]
VN <- cf$table[2,2]
table_tmp <- rbind(VP,FN)
table_tmp <- rbind(table_tmp,FP)
table_tmp <- rbind(table_tmp,VN)
table <- cbind(table, table_tmp)
write.table(overall, file="./metrics/Accuracy_Model_500bp_Class2.txt", row.names=F)
write.table(byClass, file="./metrics/byClass_Model_500bp_Class2.txt", row.names=F)
write.table(pred, file="./metrics/predictions_Model_500bp_Class2.txt", row.names=F)
write.table(table, file="./metrics/confusionMatrix_Model_500bp_Class2.txt", row.names=F)

save.image("models.RData")

