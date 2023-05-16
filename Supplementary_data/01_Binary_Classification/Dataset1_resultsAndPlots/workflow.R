library(data.table)
library(dplyr)
library(randomForest)
library(caret)
library(ranger)
library(e1071)
library(doParallel)

#load("./models/models.RData")

cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

et_grid =  expand.grid(mtry = 4:7, numRandomCuts = 1:10)
control <- trainControl(method="repeatedcv", number=10, repeats=3)
ctrl <- trainControl(method="cv", summaryFunction=twoClassSummary, classProbs=T, savePredictions = T)

#### Run classification from the best features (selected by RFE)
## Classification 1: Vector viruses and Other viruses

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/nonKmer4_13_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_genome_class1 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_genome_class1 <- train(classes ~ ., data=train, method="extraTrees", preProc=c("center", "scale"), trControl=ctrl,tuneGrid = et_grid)
#saveRDS(model_genome_class1, "model_genome_class1_roc.rds")
test <- fread(file="frag_10000bp_class1.csv", header=T, stringsAsFactors = T, sep=",")
test <- test[!duplicated(test), ]
names.use <- names(test)[(names(test) %in% col_genome_class1)]
test <- subset(test, select = c(names.use))
pred_test <- predict(model_genome_class1, newdata = test, type = "prob") %>% mutate('classes'=names(.)[apply(., 1, which.max)])
frag <- replicate(nrow(test), "10.000bp")
pred_class1 <- cbind(pred_test, frag)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/nonKmer3_12_5000bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
#subseq <- subseq[,-1]
subseq <- subset(subseq, select = -c(V1))
col_5000bp_class1 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_5000bp_class1 <- train(classes ~ ., data=train, method="extraTrees", preProc=c("center", "scale"), trControl=ctrl,tuneGrid = et_grid)
#saveRDS(model_5000bp_class1, "model_5000bp_class1_roc.rds")
test <- fread(file="frag_5000bp_class1.csv", header=T, stringsAsFactors = T, sep=",")
test <- test[!duplicated(test), ]
names.use <- names(test)[(names(test) %in% col_5000bp_class1)]
test <- subset(test, select = c(names.use))
pred_test <- predict(model_5000bp_class1, newdata = test, type = "prob") %>% mutate('classes'=names(.)[apply(., 1, which.max)])
frag <- replicate(nrow(test), "5.000bp")
pred_class1 <- cbind(pred_test, frag)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/nonKmer3_12_3000bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_3000bp_class1 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_3000bp_class1 <- train(classes ~ ., data=train, method = "knn", trControl = control, preProcess = c("center","scale"), tuneLength = 20)
#saveRDS(model_3000bp_class1, "model_3000bp_class1_roc.rds")
test <- fread(file="frag_3000bp_class1.csv", header=T, stringsAsFactors = T, sep=",")
test <- test[!duplicated(test), ]
names.use <- names(test)[(names(test) %in% col_3000bp_class1)]
test <- subset(test, select = c(names.use))
pred_test <- predict(model_3000bp_class1, newdata = test, type = "prob") %>% mutate('classes'=names(.)[apply(., 1, which.max)])
frag <- replicate(nrow(test), "3.000bp")
pred_class1 <- cbind(pred_test, frag)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/nonKmer4_13_1000bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_1000bp_class1 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_1000bp_class1 <- train(classes ~ ., data=train, method = "knn", trControl = control, preProcess = c("center","scale"), tuneLength = 20)
#saveRDS(model_1000bp_class1, "model_1000bp_class1_roc.rds")
test <- fread(file="frag_1000bp_class1.csv", header=T, stringsAsFactors = T, sep=",")
test <- test[!duplicated(test), ]
names.use <- names(test)[(names(test) %in% col_1000bp_class1)]
test <- subset(test, select = c(names.use))
pred_test <- predict(model_1000bp_class1, newdata = test, type = "prob") %>% mutate('classes'=names(.)[apply(., 1, which.max)])
frag <- replicate(nrow(test), "1.000bp")
pred_class1 <- cbind(pred_test, frag)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/nonKmer3_12_500bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_500bp_class1 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_500bp_class1 <- train(classes ~ ., data=train, method = "knn", trControl = control, preProcess = c("center","scale"), tuneLength = 20)
#saveRDS(model_500bp_class1, "model_500bp_class1_roc.rds")
test <- fread(file="frag_300bp_class1.csv", header=T, stringsAsFactors = T, sep=",")
test <- test[!duplicated(test), ]
names.use <- names(test)[(names(test) %in% col_500bp_class1)]
test <- subset(test, select = c(names.use))
pred_test <- predict(model_500bp_class1, newdata = test, type = "prob") %>% mutate('classes'=names(.)[apply(., 1, which.max)])
frag <- replicate(nrow(test), "300 bp")
pred_class1 <- cbind(pred_test, frag)

write.csv(pred_class1, file="Output_class1.csv", row.names=F) ### score file

## Classification 2: Arboviruses and Vector specific viruses

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/Kmer3_21_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_genome_class2 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_genome_class2 <- train(classes ~ ., data=train, method="extraTrees", preProc=c("center", "scale"), trControl=ctrl,tuneGrid = et_grid)
#saveRDS(model_genome_class2, "./models/model_genome_class2_roc.rds")
test <- fread(file="frag_10000bp_class2.csv", header=T, stringsAsFactors = T, sep=",")
test <- test[!duplicated(test), ]
names.use <- names(test)[(names(test) %in% col_genome_class2)]
test <- subset(test, select = c(names.use))
pred_test <- predict(model_genome_class2, newdata = test, type = "prob") %>% mutate('classes'=names(.)[apply(., 1, which.max)])
frag <- replicate(nrow(test), "10.000bp")
pred_class2 <- cbind(pred_test, frag)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/Kmer3_12_5000bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_5000bp_class2 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_5000bp_class2 <- train(classes ~ ., data=train, method="extraTrees", preProc=c("center", "scale"), trControl=ctrl,tuneGrid = et_grid)
#saveRDS(model_5000bp_class2, "model_5000bp_class2_roc.rds")
test <- fread(file="frag_5000bp_class2.csv", header=T, stringsAsFactors = T, sep=",")
test <- test[!duplicated(test), ]
names.use <- names(test)[(names(test) %in% col_5000bp_class2)]
test <- subset(test, select = c(names.use))
pred_test <- predict(model_5000bp_class2, newdata = test, type = "prob") %>% mutate('classes'=names(.)[apply(., 1, which.max)])
frag <- replicate(nrow(test), "5.000bp")
pred_class2 <- cbind(pred_test, frag)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/Kmer3_12_3000bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_3000bp_class2 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_3000bp_class2 <- train(classes ~ ., data=train, method="extraTrees", preProc=c("center", "scale"), trControl=ctrl,tuneGrid = et_grid)
#saveRDS(model_3000bp_class2, "model_3000bp_class2_roc.rds")
test <- fread(file="frag_3000bp_class2.csv", header=T, stringsAsFactors = T, sep=",")
test <- test[!duplicated(test), ]
names.use <- names(test)[(names(test) %in% col_3000bp_class2)]
test <- subset(test, select = c(names.use))
pred_test <- predict(model_3000bp_class2, newdata = test, type = "prob") %>% mutate('classes'=names(.)[apply(., 1, which.max)])
frag <- replicate(nrow(test), "3.000bp")
pred_class2 <- cbind(pred_test, frag)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/Kmer3_12_1000bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
subseq <- subset(subseq, select = -c(V1))
col_1000bp_class2 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_1000bp_class2 <- train(classes ~ ., data=train, method="extraTrees", preProc=c("center", "scale"), trControl=ctrl,tuneGrid = et_grid)
##saveRDS(model_1000bp_class2, "model_1000bp_class2_roc.rds")
test <- fread(file="frag_1000bp_class2.csv", header=T, stringsAsFactors = T, sep=",")
test <- test[!duplicated(test), ]
names.use <- names(test)[(names(test) %in% col_1000bp_class2)]
test <- subset(test, select = c(names.use))
pred_test <- predict(model_1000bp_class2, newdata = test, type = "prob") %>% mutate('classes'=names(.)[apply(., 1, which.max)])
frag <- replicate(nrow(test), "1.000bp")
pred_class2 <- cbind(pred_test, frag)

subseq <- fread(file="/disk09/uv2000/tavinbio/amanda/mnmers_recent/Fragments_models/mnmers_from_dataset/Kmer3_12_500bp_vImp.csv", header=T, stringsAsFactors = T, sep=",")
col_500bp_class2 <- colnames(subseq)
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
model_500bp_class2 <- train(classes ~ ., data=train, method = "knn", trControl = control, preProcess = c("center","scale"), tuneLength = 20)
#saveRDS(model_500bp_class2, "model_500bp_class2_roc.rds")
test <- fread(file="frag_300bp_class2.csv", header=T, stringsAsFactors = T, sep=",")
test <- test[!duplicated(test), ]
names.use <- names(test)[(names(test) %in% col_500bp_class2)]
test <- subset(test, select = c(names.use))
pred_test <- predict(model_500bp_class2, newdata = test, type = "prob") %>% mutate('classes'=names(.)[apply(., 1, which.max)])
frag <- replicate(nrow(test), "300 bp")
pred_class2 <- cbind(pred_test, frag)

write.csv(pred_class2, file="Output_class2.csv", row.names=F) ### score file

