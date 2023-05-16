library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)

sharon31 <- fread("metrics_sharon_31mer.csv", header=T)
sharon40 <- fread("metrics_sharon_40mer.csv", header=T)
sharon13 <- fread("metrics_sharon_13mer.csv", header=T)
sharon22 <- fread("metrics_sharon_22mer.csv", header=T)
sharon <- rbind(sharon31,sharon40)
sharon <- rbind(sharon,sharon13)
sharon <- rbind(sharon,sharon22)
classes <- replicate(nrow(sharon), "Sharon")
sharon <- cbind(sharon, classes)

p2 <- ggplot(sharon, aes(x=algorithm, y=adjrand, group=mn, color=mn)) +  geom_line() + geom_point() + theme_linedraw() + theme(strip.text.y.left = element_text(angle = 0))
p3 <- p2 + labs(x=NULL, y = "Adjusted rand index", color = NULL) +  theme(axis.text = element_text(colour = "black", size=11,family='Arial'))
p5 <- p3 + theme(axis.title = element_text(family = "Arial", size = 11), axis.text = element_text(family = "Arial"), axis.text.x = element_text(family = "Arial"), legend.text = element_text(size = 11, family = "Arial"))
p1 <- p3 + theme_linedraw() + theme(strip.text = element_text(colour = 'black')) +   theme(strip.background =element_rect(fill="white"))
sharon.plt <- p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values=c("{0,4}-mer"="#CC3300","{1,3}-mer"="#6BAED6","{2,2}-mer"="#4292C6","{3,1}-mer"="#2171B5"))
sharon.plt

png("sharon1_adjrand.png", width = 2000, height =1000, res = 400)
sharon.plt
dev.off()

####

mockup31 <- fread("metrics_mockup_31mer.csv", header=T)
mockup40 <- fread("metrics_mockup_40mer.csv", header=T)
mockup13 <- fread("metrics_mockup_13mer.csv", header=T)
mockup22 <- fread("metrics_mockup_22mer.csv", header=T)
mockup <- rbind(mockup31,mockup40)
mockup <- rbind(mockup,mockup13)
mockup <- rbind(mockup,mockup22)
classes <- replicate(nrow(mockup), "Mockup")
mockup <- cbind(mockup, classes)

p2 <- ggplot(mockup, aes(x=algorithm, y=adjrand, group=mn, color=mn)) +  geom_line() + geom_point() + theme_linedraw() + theme(strip.text.y.left = element_text(angle = 0))
p3 <- p2 + labs(x=NULL, y = "Adjusted rand index", color = NULL) +  theme(axis.text = element_text(colour = "black", size=11,family='Arial'))
p5 <- p3 + theme(axis.title = element_text(family = "Arial", size = 11), axis.text = element_text(family = "Arial"), axis.text.x = element_text(family = "Arial"), legend.text = element_text(size = 11, family = "Arial"))
p1 <- p3 + theme_linedraw() + theme(strip.text = element_text(colour = 'black')) +   theme(strip.background =element_rect(fill="white"))
mockup.plt <- p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values=c("{0,4}-mer"="#CC3300","{1,3}-mer"="#6BAED6","{2,2}-mer"="#4292C6","{3,1}-mer"="#2171B5"))

png("mockup_adjrand.png", width = 2000, height =1000, res = 400)
mockup.plt
dev.off()

low31 <- fread("metrics_low_31mer.csv", header=T)
low40 <- fread("metrics_low_40mer.csv", header=T)
low13 <- fread("metrics_low_13mer.csv", header=T)
low22 <- fread("metrics_low_22mer.csv", header=T)
low <- rbind(low31,low40)
low <- rbind(low,low13)
low <- rbind(low,low22)
classes <- replicate(nrow(low), "Low")
low <- cbind(low, classes)

p2 <- ggplot(low, aes(x=algorithm, y=adjrand, group=mn, color=mn)) +  geom_line() + geom_point() + theme_linedraw() + theme(strip.text.y.left = element_text(angle = 0))
p3 <- p2 + labs(x=NULL, y = "Adjusted rand index", color = NULL) +  theme(axis.text = element_text(colour = "black", size=11,family='Arial'))
p5 <- p3 + theme(axis.title = element_text(family = "Arial", size = 11), axis.text = element_text(family = "Arial"), axis.text.x = element_text(family = "Arial"), legend.text = element_text(size = 11, family = "Arial"))
p1 <- p3 + theme_linedraw() + theme(strip.text = element_text(colour = 'black')) +   theme(strip.background =element_rect(fill="white"))
low.plt <- p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values=c("{0,4}-mer"="#CC3300","{1,3}-mer"="#6BAED6","{2,2}-mer"="#4292C6","{3,1}-mer"="#2171B5"))
low.plt

png("low_adjrand.png", width = 2000, height =1000, res = 400)
low.plt
dev.off()

med31 <- fread("metrics_med_31mer.csv", header=T)
med40 <- fread("metrics_med_40mer.csv", header=T)
med13 <- fread("metrics_med_13mer.csv", header=T)
med22 <- fread("metrics_med_22mer.csv", header=T)
med <- rbind(med31,med40)
med <- rbind(med,med13)
med <- rbind(med,med22)
classes <- replicate(nrow(med), "Medium")
med <- cbind(med, classes)

p2 <- ggplot(med, aes(x=algorithm, y=adjrand, group=mn, color=mn)) +  geom_line() + geom_point() + theme_linedraw() + theme(strip.text.y.left = element_text(angle = 0))
p3 <- p2 + labs(x=NULL, y = "Adjusted rand index", color = NULL) +  theme(axis.text = element_text(colour = "black", size=11,family='Arial'))
p5 <- p3 + theme(axis.title = element_text(family = "Arial", size = 11), axis.text = element_text(family = "Arial"), axis.text.x = element_text(family = "Arial"), legend.text = element_text(size = 11, family = "Arial"))
p1 <- p3 + theme_linedraw() + theme(strip.text = element_text(colour = 'black')) +   theme(strip.background =element_rect(fill="white"))
med.plt <- p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values=c("{0,4}-mer"="#CC3300","{1,3}-mer"="#6BAED6","{2,2}-mer"="#4292C6","{3,1}-mer"="#2171B5"))
med.plt

png("med_adjrand.png", width = 2000, height =1000, res = 400)
med.plt
dev.off()

###

high31 <- fread("metrics_high_31mer.csv", header=T)
high40 <- fread("metrics_high_40mer.csv", header=T)
high13 <- fread("metrics_high_13mer.csv", header=T)
high22 <- fread("metrics_high_22mer.csv", header=T)
high <- rbind(high31,high40)
high <- rbind(high,high13)
high <- rbind(high,high22)
classes <- replicate(nrow(high), "High")
high <- cbind(high, classes)

p2 <- ggplot(high, aes(x=algorithm, y=adjrand, group=mn, color=mn)) +  geom_line() + geom_point() + theme_linedraw() + theme(strip.text.y.left = element_text(angle = 0))
p3 <- p2 + labs(x=NULL, y = "Adjusted rand index", color = NULL) +  theme(axis.text = element_text(colour = "black", size=11,family='Arial'))
p5 <- p3 + theme(axis.title = element_text(family = "Arial", size = 11), axis.text = element_text(family = "Arial"), axis.text.x = element_text(family = "Arial"), legend.text = element_text(size = 11, family = "Arial"))
p1 <- p3 + theme_linedraw() + theme(strip.text = element_text(colour = 'black')) +   theme(strip.background =element_rect(fill="white"))
med.plt <- p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values=c("{0,4}-mer"="#CC3300","{1,3}-mer"="#6BAED6","{2,2}-mer"="#4292C6","{3,1}-mer"="#2171B5"))
med.plt

png("high_adjrand.png", width = 2000, height =1000, res = 400)
high.plt
dev.off()

#### Mean

high31 <- fread("metrics_high_31mer.csv", header=T)
high40 <- fread("metrics_high_40mer.csv", header=T)
high13 <- fread("metrics_high_13mer.csv", header=T)
high22 <- fread("metrics_high_22mer.csv", header=T)
high <- rbind(high31,high40)
high <- rbind(high,high13)
high <- rbind(high,high22)
classes <- replicate(nrow(high), "High")
high <- cbind(high, classes)
high <- high[,-1]
sharon <- fread("metrics_mean_sharon.csv")
sharon <- sharon[!grepl("low", sharon$db),]
classes <- replicate(nrow(sharon), "Sharon")
sharon <- cbind(sharon, classes)

mockup <- fread("metrics_mean_mockup.csv")
mockup <- mockup[!grepl("low", mockup$db),]
classes <- replicate(nrow(mockup), "Mockup")
mockup <- cbind(mockup, classes)

low <- fread("metrics_mean_low.csv")
classes <- replicate(nrow(low), "Low")
low <- cbind(low, classes)

med <- fread("metrics_mean_med.csv")
med <- med[!grepl("low", med$db),]
classes <- replicate(nrow(med), "Medium")
med <- cbind(med, classes)

df <- rbind(sharon, mockup)
df <- rbind(df, med)
df <- rbind(df, low)
df <- rbind(df, high)
#df <- df[!grepl("Spectral",df$algorithm),]

df$classes_t <- factor(df$classes, levels = c("High", "Medium", "Low","Mockup","Sharon"))

p2 <- ggplot(df, aes(x=algorithm, y=adjrand, group=mn, color=mn)) +  geom_line() + geom_point() + theme_linedraw() + theme(strip.text.y.left = element_text(angle = 0)) +  facet_grid(classes_t~.)
p3 <- p2 + labs(x=NULL, y = "Mean Adjusted Rand Index", color = NULL) +  theme(axis.text = element_text(colour = "black", size=8,family='Arial'))
p5 <- p3 + theme(axis.title = element_text(family = "Arial", size = 8), axis.text = element_text(family = "Arial"), axis.text.x = element_text(family = "Arial"), legend.text = element_text(size = 8, family = "Arial"))
p1 <- p3 + theme_linedraw() + theme(strip.text = element_text(colour = 'black')) +   theme(strip.background =element_rect(fill="white"))
plt <- p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values=c("{0,4}-mer"="#CC3300","{1,3}-mer"="#6BAED6","{2,2}-mer"="#4292C6","{3,1}-mer"="#2171B5"))
plt

teste <- ggpar(plt, font.x = c(5,"black"),font.y = c(8,"black"),font.tickslab = c(8,"black"), legend = "top")
teste

png("Figure3_Andrade_et_al.png", width = 2000, height = 2000, res = 400)
teste
dev.off()
