# load libraries
library(ggplot2)
library(data.table)
library(ggpubr)
library(data.table)
library(mnmer)
library(dplyr)
library(doParallel)
library(doParallel)
library(microbenchmark)

set.seed(12345)

# we set only one thread
cl <- makePSOCKcluster(1)
registerDoParallel(cl)

# Function that generates Feature Matrices using our mnmer R package
# Note: You should unzip all fasta files prior to using this script.

get_train <- function(m,n){
    
    g1 <- mnmer::readNumFASTA("linear.fasta", 500, TRUE, 0.1)
    g2 <- mnmer::readNumFASTA("multisegmented.fasta", 500, TRUE, 0.1)
    g3 <- mnmer::readNumFASTA("twosegments.fasta", 500, TRUE,0.1)
    g4 <- mnmer::readNumFASTA("trisegmented.fasta", 500, TRUE,0.1)
    
    class3 <- mnmer(g3, m,n)
    classes <- replicate(nrow(class3), "twosegments")
    class3 <- cbind(class3, classes)
    
    class1 <- mnmer(g1, m,n)
    classes <- replicate(nrow(class1), "linear")
    class1 <- cbind(class1, classes)
    
    class2 <- mnmer(g2, m,n)
    classes <- replicate(nrow(class2), "multisegmented")
    class2 <- cbind(class2, classes)
    
    class4 <- mnmer(g4, m,n)
    classes <- replicate(nrow(class4), "trisegmented")
    class4 <- cbind(class4, classes)
    
    mn <- rbind(class1, class2)
    mn <- rbind(mn,class3)
    mn <- rbind(mn,class4)
    mn <- mn[,-1]
    
    mn$classes <- as.factor(mn$classes)
    train_index <- createDataPartition(mn$classes, p=0.8, list=FALSE)
    train <- mn[train_index, ]
    
    return (train)
    
}

# Benchmark function

runtime <- benchmark("{0,4}-mer"= {train <- get_train(4,0)},
                     "{0,5}-mer"= {train <- get_train(5,0)},
                     "{0,3}-mer"= {train <- get_train(3,0)},
                     "{4,1}-mer"= {train <- get_train(4,1)},
                     "{3,2}-mer"= {train <- get_train(3,2)},
                     "{2,3}-mer"= {train <- get_train(2,3)},
                     "{1,4}-mer"= {train <- get_train(1,4)},
                     "{3,1}-mer"= {train <- get_train(3,1)},
                     "{2,2}-mer"= {train <- get_train(2,2)},
                     "{1,3}-mer"= {train <- get_train(1,3)},
                     "{2,1}-mer"= {train <- get_train(2,1)},
                     "{1,2}-mer"= {train <- get_train(1,2)},
                     times = 100)
# Add class information
runtime$classes <- c("k = 3", "k = 4", "k = 5", "k = 3", "k = 4", "k = 5", "k = 3", "k = 4", "k = 5", "k = 4", "k = 5", "k = 5")

# To convert time from microseconds to seconds
runtime$elapsed <- as.numeric(as.character(runtime$elapsed)) / 1000

# Plotting 
p2 <- ggplot(runtime, aes(x=test, y=elapsed, group=classes, color=classes))  + geom_point(size=3) + geom_line(linetype="dashed", size=1) + theme_linedraw() + theme(strip.text.y.left = element_text(angle = 0))
p3 <- p2 + labs(x=NULL, y = "Elapsed time (sec)",color = NULL) +  theme(axis.text = element_text(colour = "black", size=11,family='Arial'))
p5 <- p3 + theme(axis.title = element_text(family = "Arial", size = 11), axis.text = element_text(family = "Arial"), axis.text.x = element_text(family = "Arial"), legend.text = element_text(size = 11, family = "Arial"))
p1 <- p3 + theme_linedraw() + theme(strip.text = element_text(colour = 'black')) +   theme(strip.background =element_rect(fill="white"))
plt <- p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values=c("k = 3"="#1B9E77","k = 4"="#D95F02","k = 5"="#7570B3"))
plt

plot <- ggpar(plt, font.x = c(5,"black"),font.y = c(8,"black"),font.tickslab = c(8,"black"), legend = "right")

png("/home/amanda/Supplementary_mnmer/04_Benchmarking_feature_gen/SupplementaryFigure4.png", width = 2000, height = 1000, res = 600)
plot
dev.off()

save.image("rbenchm_getmnmers.RData")