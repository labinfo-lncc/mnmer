![logo](https://user-images.githubusercontent.com/57667417/191082345-57fed066-37e9-4a8a-a65a-c9562d0625a4.png)

# Conditional frequency distribution

The (m,n)-mer R package was created to summarize biological data into numerical characteristics. It reads a FASTA file and generates a table describing the conditional frequency distribution of the selected (m,n)-mer in the sequences. This output is combined with class information to generate a feature matrix for classification. 

(m,n)-mers are an alternative case of k-mers (Figure 1). We proposed the replacement of the unconditional k-mer frequency 
by the conditional frequency, which represents the relative frequency of the n-mer conditioned to the relative frequency of m-mer plus n-mer (for more details and performance comparisson please see Andrade et al., 2022 (in press)). 

![Fig1](https://user-images.githubusercontent.com/57667417/191081859-0b0ae464-f257-4c82-9dea-8d4629605357.png)

**Fig 1.** Comparing k-mer to mn-mer relative frequency.

According to Figure 2, the k-mers are represented as (0,k) and the mn-mers as (m,n).

![Fig2](https://user-images.githubusercontent.com/57667417/191081936-1aed5ca6-9c88-4d4d-a46b-e1ccae0bcafe.png)

**Fig 2.** Numeric representation.

The output table (Figure 3) includes the fasta file accession numbers as an ID column, the relative frequency of mn-mers up to 4k columns, and class information. 

![Fig3](https://user-images.githubusercontent.com/57667417/191082016-b6835c4c-c115-498d-a2d1-c7d93ec20fe5.png)
**Fig 3.** Output example.


## Dependencies



## Installation

```
library(devtools)

devtools::install_github("labinfo-lncc/mnmer", ref="main")
```



## Quick Start: Running (m,n)-mer on example dataset

```
library("mnmer")
```

Assume we need to distinguish between viruses detected in mosquito samples and viruses that exclusively infect plants. The mn function generates the feature matrix using conditional probability from the datasets mosquito vir.fasta and plant vir.fasta. This function can generate both k-mers and mn-mers.

#### Producing k-mers

#### Producing (m,n)-mers 

The outputs are written into mosquito.csv and plant.csv files.

For classification outside of the mnmer program, we utilize the feature matrices. Here's a real-world example of code:

```

library(data.table)
library(caret)

ctrl <- trainControl(method="repeatedcv", number=10, repeats=3)
mos <- fread(file="mosquito.csv", header=T, stringsAsFactors = T, sep=",")
classes <- replicate(nrow(mos), "mosquito.vir")
r2 <- cbind(mos,classes)
plant <- fread(file="plant.csv", header=T, stringsAsFactors = T, sep=",")
classes <- replicate(nrow(plant), "plant.vir")
r1 <- cbind(plant,classes)

subseq <- rbind(r1,r2)
subseq <- subset(subseq, select = -c(seqid))
train_index <- createDataPartition(subseq$classes, p=0.8, list=FALSE)
train <- subseq[train_index, ]
test <- subseq[-train_index, ]
control <- trainControl(method="cv", summaryFunction=twoClassSummary, classProbs=T, savePredictions = T)
roc <- train(classes ~ ., data=train, method="rf", preProc=c("center"), trControl=control)

```

This classification produces the ROC curve and metrics shown below:


![Rplot01](https://user-images.githubusercontent.com/57667417/191288837-2f13cee0-96f8-48fb-a4e0-e7e28d832efe.png)


Metrics | Value
--- | ---
AUC | 1.00000
ROC | 0.99969
Sensibility | 0.97000
Specificity | 0.9980

## Citation

All data used in this project are stored at: (link) (FASTA files), (link) (Feature Matrices), and (link) (Results Metrics). 



## Work in progress

If you have any queries or find a bug, please submit an issue on GitHub or email atrv@lncc.br.
