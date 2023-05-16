library(data.table)
options(scipen=999)

### Creating an combine metrics table, using mean values

temp = list.files(pattern="byClass*")
metrics <- data.frame(matrix(ncol = 10, nrow = 1))
colnames(metrics) <- c("mn","db","algorithm","Sensitivity","Specificity","Precision","F1","Recall","Accuracy","AUC")

# Assuming your dataframe is called 'df'
num_columns <- 550
chunk_size <- 11

# Calculate the number of chunks
num_chunks <- 50

for (i in 1:length(temp)){
    byClass <-fread(paste0(temp[i]), sep = ",")
    byClass <- byClass[,-1]
    algorithm = unlist(strsplit(temp[i], split='_', fixed=TRUE))[2]
    db = unlist(strsplit(temp[i], split='_', fixed=TRUE))[3]
    mn = unlist(strsplit(temp[i], split='_', fixed=TRUE))[4]
    mn = unlist(strsplit(mn, split='.', fixed=TRUE))[1]
    chunks <- vector("list", num_chunks)
    for (i in 1:num_chunks) {
        start_col <- (i - 1) * chunk_size + 1
        end_col <- min(start_col + chunk_size - 1, num_columns)
        chunks[[i]] <- byClass[, start_col:end_col]
    }
    combined_df <- do.call(rbind, chunks)
    col_means <- colMeans(combined_df)
    col_means <- t(col_means)
    Accuracy <- fread(paste0("Accuracy_",algorithm,"_",db,"_",mn,".txt"))
    Accuracy <- Accuracy[,-1]
    Accuracy <- Accuracy[1,]
    Accuracy <- rowMeans(Accuracy)
    AUC <- fread(paste0("AUC_",db,"_",mn,".txt"))
    AUC <- AUC[-1,]
    AUC <- subset(AUC, select=c(paste0("AUC_",algorithm)))
    AUC <- colMeans(AUC)
    tmp_df <- data.frame(mn,algorithm,db,col_means,Accuracy,AUC)
    tmp_df <- subset(tmp_df, select = c(mn,db,algorithm,Sensitivity,Specificity,Precision,F1,Recall,Accuracy,AUC))
    metrics <- rbind(metrics,tmp_df)
}


metrics <- metrics[-1,]
rownames(metrics) <- NULL
write.csv(metrics, "SupplementaryTable_2.csv", row.names = F)
db4 <- metrics[grepl("dataset4", metrics$db),]
write.csv(db4, "SupplementaryTable_2_db4.csv", row.names = F)
db5 <- metrics[grepl("dataset5", metrics$db),]
write.csv(db5, "SupplementaryTable_2_db5.csv", row.names = F)

###


### Creating an combine metrics table, using mean values
setwd("/home/amanda/Downloads/03_Metrics/")
temp = list.files(pattern="byClass*")
metrics <- data.frame(matrix(ncol = 10, nrow = 1))
colnames(metrics) <- c("mn","db","algorithm","Sensitivity","Specificity","Precision","F1","Recall","Accuracy","AUC")

for (i in 1:length(temp)){
    byClass <-fread(paste0(temp[i]), sep = ",")
    byClass <- byClass[,-1]
    algorithm = unlist(strsplit(temp[i], split='_', fixed=TRUE))[2]
    db = unlist(strsplit(temp[i], split='_', fixed=TRUE))[3]
    mn = unlist(strsplit(temp[i], split='_', fixed=TRUE))[4]
    mn = unlist(strsplit(mn, split='.', fixed=TRUE))[1]
    col_means <- colMeans(byClass)
    col_means <- t(col_means)
    Accuracy <- fread(paste0("Accuracy_",algorithm,"_",db,"_",mn,".txt"))
    Accuracy <- Accuracy[,-1]
    Accuracy <- Accuracy[1,]
    Accuracy <- rowMeans(Accuracy)
    AUC <- fread(paste0("AUC_",db,"_",mn,".txt"))
    AUC <- subset(AUC, select=c(paste0("AUC_",algorithm)))
    AUC <- colMeans(AUC)
    tmp_df <- data.frame(mn,algorithm,db,col_means,Accuracy,AUC)
    tmp_df <- subset(tmp_df, select = c(mn,db,algorithm,Sensitivity,Specificity,Precision,F1,Recall,Accuracy,AUC))
    metrics <- rbind(metrics,tmp_df)
}


metrics <- metrics[-1,]
rownames(metrics) <- NULL
write.csv(metrics, "SupplementaryTable_2.csv", row.names = F)
db4 <- metrics[grepl("dataset4", metrics$db),]
write.csv(db4, "SupplementaryTable_2_db4.csv", row.names = F)
db5 <- metrics[grepl("dataset5", metrics$db),]
write.csv(db5, "SupplementaryTable_2_db5.csv", row.names = F)
db5 <- metrics[grepl("dataset6", metrics$db),]
write.csv(db5, "SupplementaryTable_2_db6.csv", row.names = F)
###

