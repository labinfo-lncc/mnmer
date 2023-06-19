# Load required libraries
library(data.table)
library(dplyr)
library(randomForest)
library(caret)
library(ranger)
library(e1071)
library(mnmer)
library(kknn)
library(mda)
library(HDclassif)
library(sda)
library(pROC)
library(doParallel)

# Set up parallel processing
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

# Set the seed for reproducibility
set.seed(12345)

# Set the train control method for cross-validation
control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = multiClassSummary)
metric <- "AUC"

# Initialize variables for storing results. These variables will store the results of each iteration for different evaluation metrics and models.
overall_knn <- c(1:7)
overall_knn <- as.data.frame(overall_knn)

overall_pda <- c(1:7)
overall_pda <- as.data.frame(overall_pda)

overall_sda <- c(1:7)
overall_sda <- as.data.frame(overall_sda)

overall_rf <- c(1:7)
overall_rf <- as.data.frame(overall_rf)

overall_C5 <- c(1:7)
overall_C5 <- as.data.frame(overall_C5)

byClass_knn <- c(1:4)
byClass_knn <- as.data.frame(byClass_knn)

byClass_pda <- c(1:4)
byClass_pda <- as.data.frame(byClass_pda)

byClass_sda <- c(1:4)
byClass_sda <- as.data.frame(byClass_sda)

byClass_rf <- c(1:4)
byClass_rf <- as.data.frame(byClass_rf)

byClass_C5 <- c(1:4)
byClass_C5 <- as.data.frame(byClass_C5)

AUC_knn <- c(1:1)
AUC_knn <- as.data.frame(AUC_knn)

AUC_pda <- c(1:1)
AUC_pda <- as.data.frame(AUC_pda)

AUC_sda <- c(1:1)
AUC_sda <- as.data.frame(AUC_sda)

AUC_rf <- c(1:1)
AUC_rf <- as.data.frame(AUC_rf)

AUC_C5 <- c(1:1)
AUC_C5 <- as.data.frame(AUC_C5)

# Perform iterations
for (i in 1:50) {
  # Read random FASTA files
  # g1, g2, g3, and g4 are different sets of sequences
  # read from FASTA files
  
  # Generate feature matrix for g3
  # class3 will contain the generated features and assigned class labels
  class3 <- mnmer(g3, 3, 1)
  classes <- replicate(nrow(class3), "Amarillovirales")
  class3 <- cbind(class3, classes)

  # Generate feature matrix for g1
  class1 <- mnmer(g1, 3, 1)
  classes <- replicate(nrow(class1), "Bunyavirales")
  class1 <- cbind(class1, classes)

  # Generate feature matrix for g2
  class2 <- mnmer(g2, 3, 1)
  classes <- replicate(nrow(class2), "Nidovirales")
  class2 <- cbind(class2, classes)

  # Generate feature matrix for g4
  class4 <- mnmer(g4, 3, 1)
  classes <- replicate(nrow(class4), "Mononegavirales")
  class4 <- cbind(class4, classes)

# Combine all feature matrices into a single matrix
  mn <- rbind(class1, class2)
  mn <- rbind(mn, class3)
  mn <- rbind(mn, class4)
  mn <- mn[,-1]

#### split into train and test sets

# Convert the 'classes' column to a factor
  mn$classes <- as.factor(mn$classes)

# Create an index for splitting the data into training and test sets
  train_index <- createDataPartition(mn$classes, p=0.8, list=FALSE)

# Split the data into training and test sets based on the index
  train <- mn[train_index, ]
  test <- mn[-train_index, ]

# Convert the 'classes' column in the training set to a factor
  train$classes <- as.factor(train$classes)

#### models

# Train a k-Nearest Neighbors model for 31mer
  mn_knn31 <- train(classes ~ ., data=train, method="kknn", metric=metric, trControl=control, preProcess=c("center", "scale"))
  p1 <- predict(mn_knn31, test, preProcess=c("center", "scale"), type="prob")
  AUC_tmp <- auc(multiclass.roc(test$classes, p1))
  AUC_knn <- rbind(AUC_knn, AUC_tmp)
  p1 <- predict(mn_knn31, newdata=test)
  cf <- confusionMatrix(p1, as.factor(test$classes))
  byClass_knn_tmp <- cf$byClass
  overall_knn_tmp <- cf$overall
  byClass_knn <- cbind(byClass_knn, byClass_knn_tmp)
  overall_knn <- cbind(overall_knn, overall_knn_tmp)

# Train a Probabilistic Discriminant Analysis model
  mn_pda31 <- train(classes ~ ., data=train, method="pda", metric=metric, trControl=control, preProcess=c("center", "scale"))
  p2 <- predict(mn_pda31, test, preProcess=c("center", "scale"), type="prob")
  AUC_tmp <- auc(multiclass.roc(test$classes, p2))
  AUC_pda <- rbind(AUC_pda, AUC_tmp)
  p2 <- predict(mn_pda31, newdata=test)
  cf <- confusionMatrix(p2, as.factor(test$classes))
  byClass_pda_tmp <- cf$byClass
  overall_pda_tmp <- cf$overall
  byClass_pda <- cbind(byClass_pda, byClass_pda_tmp)
  overall_pda <- cbind(overall_pda, overall_pda_tmp)

# Train a Shrinkage Discriminant Analysis model
  mn_sda31 <- train(classes ~ ., data=train, method="sda", metric=metric, trControl=control, preProcess=c("center", "scale"))
  p3 <- predict(mn_sda31, test, preProcess=c("center", "scale"), type="prob")
  AUC_tmp <- auc(multiclass.roc(test$classes, p3))
  AUC_sda <- rbind(AUC_sda, AUC_tmp)
  p3 <- predict(mn_sda31, newdata=test)
  cf <- confusionMatrix(p3, as.factor(test$classes))
  byClass_sda_tmp <- cf$byClass
  overall_sda_tmp <- cf$overall
  byClass_sda <- cbind(byClass_sda, byClass_sda_tmp)
  overall_sda <- cbind(overall_sda, overall_sda_tmp)

# Train a Random Forest model
  mn_rf31 <- train(classes ~ ., data=train, method="rf", metric=metric, trControl=control, preProcess=c("center", "scale"))
  p4 <- predict(mn_rf31, test, preProcess=c("center", "scale"), type="prob")
  AUC_tmp <- auc(multiclass.roc(test$classes, p4))
  AUC_rf <- rbind(AUC_rf, AUC_tmp)
  p4 <- predict(mn_rf31, newdata=test)
  cf <- confusionMatrix(p4, as.factor(test$classes))
  byClass_rf_tmp <- cf$byClass
  overall_rf_tmp <- cf$overall
  byClass_rf <- cbind(byClass_rf, byClass_rf_tmp)
  overall_rf <- cbind(overall_rf, overall_rf_tmp)

# Train a C5.0 Decision Tree model
  mn_C531 <- train(classes ~ ., data=train, method="C5.0Tree", metric=metric, trControl=control, preProcess=c("center", "scale"))
  p5 <- predict(mn_C531, test, preProcess=c("center", "scale"), type="prob")
  AUC_tmp <- auc(multiclass.roc(test$classes, p5))
  AUC_C5 <- rbind(AUC_C5, AUC_tmp)
  p5 <- predict(mn_C531, newdata=test)
  cf <- confusionMatrix(p5, as.factor(test$classes))
  byClass_C5_tmp <- cf$byClass
  overall_C5_tmp <- cf$overall
  byClass_C5 <- cbind(byClass_C5, byClass_C5_tmp)
  overall_C5 <- cbind(overall_C5, overall_C5_tmp)
	
  # Remove unnecessary variables
  rm(train,test, class1, class2, arbo, most)
	
  # Print progress
  print(paste0("OK run31: ", i))

}

# Combine the AUC values for different models

AUC <- cbind(AUC_knn, AUC_pda)
AUC <- cbind(AUC_rf, AUC)
AUC <- cbind(AUC_sda, AUC)
AUC <- cbind(AUC_C5, AUC)

### write tables
# Write the results to CSV files

write.csv(AUC, file="../03_Metrics/AUC_dataset4_31.txt", row.names=F)
write.csv(overall_knn, file="../03_Metrics/Accuracy_knn_dataset4_31.txt", row.names=F)
write.csv(overall_sda, file="../03_Metrics/Accuracy_sda_dataset4_31.txt", row.names=F)
write.csv(overall_pda, file="../03_Metrics/Accuracy_pda_dataset4_31.txt", row.names=F)
write.csv(overall_rf, file="../03_Metrics/Accuracy_rf_dataset4_31.txt", row.names=F)
write.csv(overall_C5, file="../03_Metrics/Accuracy_C5_dataset4_31.txt", row.names=F)
write.csv(byClass_knn, file="../03_Metrics/byClass_knn_dataset4_31.txt", row.names=F)
write.csv(byClass_sda, file="../03_Metrics/byClass_sda_dataset4_31.txt", row.names=F)
write.csv(byClass_pda, file="../03_Metrics/byClass_pda_dataset4_31.txt", row.names=F)
write.csv(byClass_rf, file="../03_Metrics/byClass_rf_dataset4_31.txt", row.names=F)
write.csv(byClass_C5, file="../03_Metrics/byClass_C5_dataset4_31.txt", row.names=F)

# Remove unnecessary variables
rm(AUC, AUC_knn, AUC_pda, AUC_rf, AUC_sda, AUC_C5, class1, class2)

#save.image("dataset4.RData")

### init the variables

overall_knn <-c(1:7)
overall_knn <- as.data.frame(overall_knn)

overall_pda <-c(1:7)
overall_pda <- as.data.frame(overall_pda)

overall_sda <-c(1:7)
overall_sda <- as.data.frame(overall_sda)

overall_rf <-c(1:7)
overall_rf <- as.data.frame(overall_rf)

overall_C5 <-c(1:7)
overall_C5 <- as.data.frame(overall_C5)

byClass_knn <- c(1:4)
byClass_knn <- as.data.frame(byClass_knn)

byClass_pda <- c(1:4)
byClass_pda <- as.data.frame(byClass_pda)

byClass_sda <- c(1:4)
byClass_sda <- as.data.frame(byClass_sda)

byClass_rf <- c(1:4)
byClass_rf <- as.data.frame(byClass_rf)

byClass_C5 <- c(1:4)
byClass_C5 <- as.data.frame(byClass_C5)

AUC_knn <-c(1:1)
AUC_knn <- as.data.frame(AUC_knn)

AUC_pda <-c(1:1)
AUC_pda <- as.data.frame(AUC_pda)

AUC_sda <-c(1:1)
AUC_sda <- as.data.frame(AUC_sda)

AUC_rf <-c(1:1)
AUC_rf <- as.data.frame(AUC_rf)

AUC_C5 <-c(1:1)
AUC_C5 <- as.data.frame(AUC_C5)

#####

for (i in 1:50){

#### read random FASTA

        g1 <- mnmer::readNumFASTA("../01_Raw_data/dataset4/Bunyavirales.fasta",500,TRUE,0.02)
        g2 <- mnmer::readNumFASTA("../01_Raw_data/dataset4/Nidovirales.fasta",500,TRUE,0.02)
        g3 <- mnmer::readNumFASTA("../01_Raw_data/dataset4/Amarillovirales.fasta",500,TRUE,0.02)
        g4 <- mnmer::readNumFASTA("../01_Raw_data/dataset4/Mononegavirales.fasta",500,TRUE,0.02)
	
#### feature matrix
	
	class3 <- mnmer(g3, 4,0)
	classes <- replicate(nrow(class3), "Amarillovirales")
	class3 <- cbind(class3, classes)

	class1 <- mnmer(g1, 4,0)
	classes <- replicate(nrow(class1), "Bunyavirales")
	class1 <- cbind(class1, classes)

	class2 <- mnmer(g2, 4,0)
	classes <- replicate(nrow(class2), "Nidovirales")
	class2 <- cbind(class2, classes)
	
	class4 <- mnmer(g4, 4,0)
	classes <- replicate(nrow(class4), "Mononegavirales")
	class4 <- cbind(class4, classes)
	
	mn <- rbind(class1, class2)
	mn <- rbind(mn,class3)
	mn <- rbind(mn,class4)
	mn <- mn[,-1]
	
#### split into train and test sets
	
	mn$classes <- as.factor(mn$classes)
	train_index <- createDataPartition(mn$classes, p=0.8, list=FALSE)
	train <- mn[train_index, ]
	test <- mn[-train_index, ]
	train$classes <- as.factor(train$classes)
	
#### models

	mn_knn40 <- train(classes~., data=train, method="kknn", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p1 <- predict(mn_knn40, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p1))
	AUC_knn <- rbind(AUC_knn, AUC_tmp)
	p1 <- predict(mn_knn40, newdata = test)
	cf <- confusionMatrix(p1, as.factor(test$classes))
	byClass_knn_tmp <- cf$byClass
	overall_knn_tmp <- cf$overall
	byClass_knn <- cbind(byClass_knn, byClass_knn_tmp)
	overall_knn <- cbind(overall_knn, overall_knn_tmp)

	mn_pda40 <- train(classes~., data=train, method="pda", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p2 <- predict(mn_pda40, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p2))
	AUC_pda <- rbind(AUC_pda, AUC_tmp)
	p2 <- predict(mn_pda40, newdata = test)
	cf <- confusionMatrix(p2, as.factor(test$classes))
	byClass_pda_tmp <- cf$byClass
	overall_pda_tmp <- cf$overall
	byClass_pda <- cbind(byClass_pda, byClass_pda_tmp)
	overall_pda <- cbind(overall_pda, overall_pda_tmp)

	mn_sda40 <- train(classes~., data=train, method="sda", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p3 <- predict(mn_sda40, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p3))
	AUC_sda <- rbind(AUC_sda, AUC_tmp)
	p3 <- predict(mn_sda40, newdata = test)
	cf <- confusionMatrix(p3, as.factor(test$classes))
	byClass_sda_tmp <- cf$byClass
	overall_sda_tmp <- cf$overall
	byClass_sda <- cbind(byClass_sda, byClass_sda_tmp)
	overall_sda <- cbind(overall_sda, overall_sda_tmp)

	mn_rf40 <- train(classes~., data=train, method="rf", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p4 <- predict(mn_rf40, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p4))
	AUC_rf <- rbind(AUC_rf, AUC_tmp)
	p4 <- predict(mn_rf40, newdata = test)
	cf <- confusionMatrix(p4, as.factor(test$classes))
	byClass_rf_tmp <- cf$byClass
	overall_rf_tmp <- cf$overall
	byClass_rf <- cbind(byClass_rf, byClass_rf_tmp)
	overall_rf <- cbind(overall_rf, overall_rf_tmp)

	mn_C540 <- train(classes~., data=train, method="C5.0Tree", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p5 <- predict(mn_C540, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p5))
	AUC_C5 <- rbind(AUC_C5, AUC_tmp)
	p5 <- predict(mn_C540, newdata = test)
	cf <- confusionMatrix(p5, as.factor(test$classes))
	byClass_C5_tmp <- cf$byClass
	overall_C5_tmp <- cf$overall
	byClass_C5 <- cbind(byClass_C5, byClass_C5_tmp)
	overall_C5 <- cbind(overall_C5, overall_C5_tmp)


	rm(train,test, class1, class2, arbo, mos)
	print(paste0("OK run40: ", i))

}

#### bind results

AUC <- cbind(AUC_knn, AUC_pda)
AUC <- cbind(AUC_rf, AUC)
AUC <- cbind(AUC_sda, AUC)
AUC <- cbind(AUC_C5, AUC)

### write tables
write.csv(AUC, file="../03_Metrics/AUC_dataset4_40.txt", row.names=F)
write.csv(overall_knn, file="../03_Metrics/Accuracy_knn_dataset4_40.txt", row.names=F)
write.csv(overall_sda, file="../03_Metrics/Accuracy_sda_dataset4_40.txt", row.names=F)
write.csv(overall_pda, file="../03_Metrics/Accuracy_pda_dataset4_40.txt", row.names=F)
write.csv(overall_rf, file="../03_Metrics/Accuracy_rf_dataset4_40.txt", row.names=F)
write.csv(overall_C5, file="../03_Metrics/Accuracy_C5_dataset4_40.txt", row.names=F)
write.csv(byClass_knn, file="../03_Metrics/byClass_knn_dataset4_40.txt", row.names=F)
write.csv(byClass_sda, file="../03_Metrics/byClass_sda_dataset4_40.txt", row.names=F)
write.csv(byClass_pda, file="../03_Metrics/byClass_pda_dataset4_40.txt", row.names=F)
write.csv(byClass_rf, file="../03_Metrics/byClass_rf_dataset4_40.txt", row.names=F)
write.csv(byClass_C5, file="../03_Metrics/byClass_C5_dataset4_40.txt", row.names=F)

rm(AUC, AUC_knn, AUC_pda, AUC_rf, AUC_sda, AUC_C5, class1, class2)

#save.image("dataset4.RData")

### init the variables

overall_knn <-c(1:7)
overall_knn <- as.data.frame(overall_knn)

overall_pda <-c(1:7)
overall_pda <- as.data.frame(overall_pda)

overall_sda <-c(1:7)
overall_sda <- as.data.frame(overall_sda)

overall_rf <-c(1:7)
overall_rf <- as.data.frame(overall_rf)

overall_C5 <-c(1:7)
overall_C5 <- as.data.frame(overall_C5)

byClass_knn <- c(1:4)
byClass_knn <- as.data.frame(byClass_knn)

byClass_pda <- c(1:4)
byClass_pda <- as.data.frame(byClass_pda)

byClass_sda <- c(1:4)
byClass_sda <- as.data.frame(byClass_sda)

byClass_rf <- c(1:4)
byClass_rf <- as.data.frame(byClass_rf)

byClass_C5 <- c(1:4)
byClass_C5 <- as.data.frame(byClass_C5)

AUC_knn <-c(1:1)
AUC_knn <- as.data.frame(AUC_knn)

AUC_pda <-c(1:1)
AUC_pda <- as.data.frame(AUC_pda)

AUC_sda <-c(1:1)
AUC_sda <- as.data.frame(AUC_sda)

AUC_rf <-c(1:1)
AUC_rf <- as.data.frame(AUC_rf)

AUC_C5 <-c(1:1)
AUC_C5 <- as.data.frame(AUC_C5)

#####

for (i in 1:50){

#### read random FASTA

        g1 <- mnmer::readNumFASTA("../01_Raw_data/dataset4/Bunyavirales.fasta",500,TRUE,0.02)
        g2 <- mnmer::readNumFASTA("../01_Raw_data/dataset4/Nidovirales.fasta",500,TRUE,0.02)
        g3 <- mnmer::readNumFASTA("../01_Raw_data/dataset4/Amarillovirales.fasta",500,TRUE,0.02)
        g4 <- mnmer::readNumFASTA("../01_Raw_data/dataset4/Mononegavirales.fasta",500,TRUE,0.02)

#### feature matrix
	
	class3 <- mnmer(g3, 2,2)
	classes <- replicate(nrow(class3), "Amarillovirales")
	class3 <- cbind(class3, classes)

	class1 <- mnmer(g1, 2,2)
	classes <- replicate(nrow(class1), "Bunyavirales")
	class1 <- cbind(class1, classes)

	class2 <- mnmer(g2, 2,2)
	classes <- replicate(nrow(class2), "Nidovirales")
	class2 <- cbind(class2, classes)
	
	class4 <- mnmer(g4, 2,2)
	classes <- replicate(nrow(class4), "Mononegavirales")
	class4 <- cbind(class4, classes)
	
	mn <- rbind(class1, class2)
	mn <- rbind(mn,class3)
	mn <- rbind(mn,class4)
	mn <- mn[,-1]
	
#### split into train and test sets
	
	mn$classes <- as.factor(mn$classes)
	train_index <- createDataPartition(mn$classes, p=0.8, list=FALSE)
	train <- mn[train_index, ]
	test <- mn[-train_index, ]
	train$classes <- as.factor(train$classes)
	
#### models

	mn_knn22 <- train(classes~., data=train, method="kknn", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p1 <- predict(mn_knn22, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p1))
	AUC_knn <- rbind(AUC_knn, AUC_tmp)
	p1 <- predict(mn_knn22, newdata = test)
	cf <- confusionMatrix(p1, as.factor(test$classes))
	byClass_knn_tmp <- cf$byClass
	overall_knn_tmp <- cf$overall
	byClass_knn <- cbind(byClass_knn, byClass_knn_tmp)
	overall_knn <- cbind(overall_knn, overall_knn_tmp)

	mn_pda22 <- train(classes~., data=train, method="pda", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p2 <- predict(mn_pda22, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p2))
	AUC_pda <- rbind(AUC_pda, AUC_tmp)
	p2 <- predict(mn_pda22, newdata = test)
	cf <- confusionMatrix(p2, as.factor(test$classes))
	byClass_pda_tmp <- cf$byClass
	overall_pda_tmp <- cf$overall
	byClass_pda <- cbind(byClass_pda, byClass_pda_tmp)
	overall_pda <- cbind(overall_pda, overall_pda_tmp)

	mn_sda22 <- train(classes~., data=train, method="sda", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p3 <- predict(mn_sda22, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p3))
	AUC_sda <- rbind(AUC_sda, AUC_tmp)
	p3 <- predict(mn_sda22, newdata = test)
	cf <- confusionMatrix(p3, as.factor(test$classes))
	byClass_sda_tmp <- cf$byClass
	overall_sda_tmp <- cf$overall
	byClass_sda <- cbind(byClass_sda, byClass_sda_tmp)
	overall_sda <- cbind(overall_sda, overall_sda_tmp)

	mn_rf22 <- train(classes~., data=train, method="rf", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p4 <- predict(mn_rf22, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p4))
	AUC_rf <- rbind(AUC_rf, AUC_tmp)
	p4 <- predict(mn_rf22, newdata = test)
	cf <- confusionMatrix(p4, as.factor(test$classes))
	byClass_rf_tmp <- cf$byClass
	overall_rf_tmp <- cf$overall
	byClass_rf <- cbind(byClass_rf, byClass_rf_tmp)
	overall_rf <- cbind(overall_rf, overall_rf_tmp)

	mn_C522 <- train(classes~., data=train, method="C5.0Tree", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p5 <- predict(mn_C522, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p5))
	AUC_C5 <- rbind(AUC_C5, AUC_tmp)
	p5 <- predict(mn_C522, newdata = test)
	cf <- confusionMatrix(p5, as.factor(test$classes))
	byClass_C5_tmp <- cf$byClass
	overall_C5_tmp <- cf$overall
	byClass_C5 <- cbind(byClass_C5, byClass_C5_tmp)
	overall_C5 <- cbind(overall_C5, overall_C5_tmp)


	rm(train,test, class1, class2, arbo, mos)
	print(paste0("OK run22: ", i))

}

#### bind results

AUC <- cbind(AUC_knn, AUC_pda)
AUC <- cbind(AUC_rf, AUC)
AUC <- cbind(AUC_sda, AUC)
AUC <- cbind(AUC_C5, AUC)

### write tables
write.csv(AUC, file="../03_Metrics/AUC_dataset4_22.txt", row.names=F)
write.csv(overall_knn, file="../03_Metrics/Accuracy_knn_dataset4_22.txt", row.names=F)
write.csv(overall_sda, file="../03_Metrics/Accuracy_sda_dataset4_22.txt", row.names=F)
write.csv(overall_pda, file="../03_Metrics/Accuracy_pda_dataset4_22.txt", row.names=F)
write.csv(overall_rf, file="../03_Metrics/Accuracy_rf_dataset4_22.txt", row.names=F)
write.csv(overall_C5, file="../03_Metrics/Accuracy_C5_dataset4_22.txt", row.names=F)
write.csv(byClass_knn, file="../03_Metrics/byClass_knn_dataset4_22.txt", row.names=F)
write.csv(byClass_sda, file="../03_Metrics/byClass_sda_dataset4_22.txt", row.names=F)
write.csv(byClass_pda, file="../03_Metrics/byClass_pda_dataset4_22.txt", row.names=F)
write.csv(byClass_rf, file="../03_Metrics/byClass_rf_dataset4_22.txt", row.names=F)
write.csv(byClass_C5, file="../03_Metrics/byClass_C5_dataset4_22.txt", row.names=F)
rm(AUC, AUC_knn, AUC_pda, AUC_rf, AUC_sda, AUC_C5, class1, class2)

#save.image("dataset4.RData")

### init the variables

overall_knn <-c(1:7)
overall_knn <- as.data.frame(overall_knn)

overall_pda <-c(1:7)
overall_pda <- as.data.frame(overall_pda)

overall_sda <-c(1:7)
overall_sda <- as.data.frame(overall_sda)

overall_rf <-c(1:7)
overall_rf <- as.data.frame(overall_rf)

overall_C5 <-c(1:7)
overall_C5 <- as.data.frame(overall_C5)

byClass_knn <- c(1:4)
byClass_knn <- as.data.frame(byClass_knn)

byClass_pda <- c(1:4)
byClass_pda <- as.data.frame(byClass_pda)

byClass_sda <- c(1:4)
byClass_sda <- as.data.frame(byClass_sda)

byClass_rf <- c(1:4)
byClass_rf <- as.data.frame(byClass_rf)

byClass_C5 <- c(1:4)
byClass_C5 <- as.data.frame(byClass_C5)

AUC_knn <-c(1:1)
AUC_knn <- as.data.frame(AUC_knn)

AUC_pda <-c(1:1)
AUC_pda <- as.data.frame(AUC_pda)

AUC_sda <-c(1:1)
AUC_sda <- as.data.frame(AUC_sda)

AUC_rf <-c(1:1)
AUC_rf <- as.data.frame(AUC_rf)

AUC_C5 <-c(1:1)
AUC_C5 <- as.data.frame(AUC_C5)

#####

for (i in 1:50){

#### read random FASTA

	g1 <- mnmer::readNumFASTA("../01_Raw_data/dataset4/Bunyavirales.fasta",500,TRUE,0.02)
	g2 <- mnmer::readNumFASTA("../01_Raw_data/dataset4/Nidovirales.fasta",500,TRUE,0.02)
	g3 <- mnmer::readNumFASTA("../01_Raw_data/dataset4/Amarillovirales.fasta",500,TRUE,0.02)
	g4 <- mnmer::readNumFASTA("../01_Raw_data/dataset4/Mononegavirales.fasta",500,TRUE,0.02)

# feature matrix
	
	class3 <- mnmer(g3, 1,3)
	classes <- replicate(nrow(class3), "Amarillovirales")
	class3 <- cbind(class3, classes)

	class1 <- mnmer(g1, 1,3)
	classes <- replicate(nrow(class1), "Bunyavirales")
	class1 <- cbind(class1, classes)

	class2 <- mnmer(g2, 1,3)
	classes <- replicate(nrow(class2), "Nidovirales")
	class2 <- cbind(class2, classes)
	
	class4 <- mnmer(g4, 1,3)
	classes <- replicate(nrow(class4), "Mononegavirales")
	class4 <- cbind(class4, classes)
	
	mn <- rbind(class1, class2)
	mn <- rbind(mn,class3)
	mn <- rbind(mn,class4)
	mn <- mn[,-1]
	
#### split into train and test sets
	
	mn$classes <- as.factor(mn$classes)
	train_index <- createDataPartition(mn$classes, p=0.8, list=FALSE)
	train <- mn[train_index, ]
	test <- mn[-train_index, ]
	train$classes <- as.factor(train$classes)
	
#### models

	mn_knn13 <- train(classes~., data=train, method="kknn", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p1 <- predict(mn_knn13, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p1))
	AUC_knn <- rbind(AUC_knn, AUC_tmp)
	p1 <- predict(mn_knn13, newdata = test)
	cf <- confusionMatrix(p1, as.factor(test$classes))
	byClass_knn_tmp <- cf$byClass
	overall_knn_tmp <- cf$overall
	byClass_knn <- cbind(byClass_knn, byClass_knn_tmp)
	overall_knn <- cbind(overall_knn, overall_knn_tmp)

	mn_pda13 <- train(classes~., data=train, method="pda", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p2 <- predict(mn_pda13, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p2))
	AUC_pda <- rbind(AUC_pda, AUC_tmp)
	p2 <- predict(mn_pda13, newdata = test)
	cf <- confusionMatrix(p2, as.factor(test$classes))
	byClass_pda_tmp <- cf$byClass
	overall_pda_tmp <- cf$overall
	byClass_pda <- cbind(byClass_pda, byClass_pda_tmp)
	overall_pda <- cbind(overall_pda, overall_pda_tmp)

	mn_sda13 <- train(classes~., data=train, method="sda", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p3 <- predict(mn_sda13, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p3))
	AUC_sda <- rbind(AUC_sda, AUC_tmp)
	p3 <- predict(mn_sda13, newdata = test)
	cf <- confusionMatrix(p3, as.factor(test$classes))
	byClass_sda_tmp <- cf$byClass
	overall_sda_tmp <- cf$overall
	byClass_sda <- cbind(byClass_sda, byClass_sda_tmp)
	overall_sda <- cbind(overall_sda, overall_sda_tmp)

	mn_rf13 <- train(classes~., data=train, method="rf", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p4 <- predict(mn_rf13, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p4))
	AUC_rf <- rbind(AUC_rf, AUC_tmp)
	p4 <- predict(mn_rf13, newdata = test)
	cf <- confusionMatrix(p4, as.factor(test$classes))
	byClass_rf_tmp <- cf$byClass
	overall_rf_tmp <- cf$overall
	byClass_rf <- cbind(byClass_rf, byClass_rf_tmp)
	overall_rf <- cbind(overall_rf, overall_rf_tmp)

	mn_C513 <- train(classes~., data=train, method="C5.0Tree", metric=metric,trControl=control, preProcess = c("center", "scale") )
	p5 <- predict(mn_C513, test, preProcess = c("center", "scale"),type = "prob")
	AUC_tmp <- auc(multiclass.roc(test$classes, p5))
	AUC_C5 <- rbind(AUC_C5, AUC_tmp)
	p5 <- predict(mn_C513, newdata = test)
	cf <- confusionMatrix(p5, as.factor(test$classes))
	byClass_C5_tmp <- cf$byClass
	overall_C5_tmp <- cf$overall
	byClass_C5 <- cbind(byClass_C5, byClass_C5_tmp)
	overall_C5 <- cbind(overall_C5, overall_C5_tmp)


	rm(train,test, class1, class2, arbo, mos)
	print(paste0("OK run13: ", i))

}

#### bind results

AUC <- cbind(AUC_knn, AUC_pda)
AUC <- cbind(AUC_rf, AUC)
AUC <- cbind(AUC_sda, AUC)
AUC <- cbind(AUC_C5, AUC)

### write tables
write.csv(AUC, file="../03_Metrics/AUC_dataset4_13.txt", row.names=F)
write.csv(overall_knn, file="../03_Metrics/Accuracy_knn_dataset4_13.txt", row.names=F)
write.csv(overall_sda, file="../03_Metrics/Accuracy_sda_dataset4_13.txt", row.names=F)
write.csv(overall_pda, file="../03_Metrics/Accuracy_pda_dataset4_13.txt", row.names=F)
write.csv(overall_rf, file="../03_Metrics/Accuracy_rf_dataset4_13.txt", row.names=F)
write.csv(overall_C5, file="../03_Metrics/Accuracy_C5_dataset4_13.txt", row.names=F)
write.csv(byClass_knn, file="../03_Metrics/byClass_knn_dataset4_13.txt", row.names=F)
write.csv(byClass_sda, file="../03_Metrics/byClass_sda_dataset4_13.txt", row.names=F)
write.csv(byClass_pda, file="../03_Metrics/byClass_pda_dataset4_13.txt", row.names=F)
write.csv(byClass_rf, file="../03_Metrics/byClass_rf_dataset4_13.txt", row.names=F)
write.csv(byClass_C5, file="../03_Metrics/byClass_C5_dataset4_13.txt", row.names=F)

rm(AUC, AUC_knn, AUC_pda, AUC_rf, AUC_sda, AUC_C5, class1, class2)

save.image("dataset4.RData")
