library(ggpubr)
library(ggplot2)
library(data.table)

#### plot for datasetA

df <- fread(file="AUCS_datasetA.csv",header=T)
df$frag = factor(df$frag, levels=c('300 bp','1.000 bp','3.000 bp','5.000 bp','10.000 bp'))

p2 <- ggplot(df, aes(x=frag, y=mean, group=type, color=type)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, position=position_dodge(0.05)) + 
  geom_line() + geom_point() + theme_linedraw() + 
  theme(strip.text.y = element_text(angle = 0))
p2
p3 <- p2 + 
  labs(x=NULL, y = "mean AUC", fill = NULL, color=NULL) + 
  theme(axis.text = element_text(colour = "black", size=12,family='Arial')) 
p3
p5 <- p3 + 
  theme(axis.title = element_text(family = "Arial", size = 12), axis.text = element_text(family = "Arial"), axis.text.x = element_text(family = "Arial"), legend.text = element_text(size = 12, family = "Arial"))+
  scale_color_manual(values=c("{0,k}-mer"="#CC3300","{1,k-1}-mer"="#6BAED6","{2,k-2}-mer"="#4292C6","{3,k-3}-mer"="#2171B5"))
p5

p1 <- p5 + 
  theme_linedraw() + 
  theme(strip.text = element_text(colour = 'black')) + 
  theme(strip.background =element_rect(fill="white")) + theme(legend.position="top") 
p1
db1 <- ggpar(p1, font.x = c(10,"black"), font.y = c(10,"black"),font.y.tickslab = c(4,"black"),legend = "top")
db1

png("AUCS_dataset1.png", width = 2000, height = 550, res = 600)
db1
dev.off()

### plots for datasetB
df <- fread(file="AUCS_datasetB.csv",header=T)
df$frag_t = factor(df$frag, levels=c('300 bp','1.000 bp','3.000 bp','5.000 bp','10.000 bp'))

p2 <- ggplot(df, aes(x=frag_t, y=mean, group=type, color=type)) +  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, position=position_dodge(0.05)) + 
  geom_line() + geom_point() + theme_linedraw() + facet_grid(k ~ .) + theme(strip.text.y.left = element_text(angle = 0))
p3 <- p2 + labs(x=NULL, y = "mean AUC", fill = NULL) +  theme(axis.text = element_text(colour = "black", size=11,family='Arial')) 
p5 <- p3 + theme(axis.title = element_text(family = "Arial", size = 11), axis.text = element_text(family = "Arial"), axis.text.x = element_text(family = "Arial"), legend.text = element_text(size = 11, family = "Arial"))+
  scale_color_manual(values=c("{0,k}-mer"="#1B9E77","{1,k-1}-mer"="#9ECAE1","{2,k-2}-mer"="#4292C6","{3,k-3}-mer"="#084594","{4,k-4}-mer"="#DEEBF7","{5,k-5}-mer"="#C6DBEF"))
p1 <- p3 + theme_linedraw() + theme(strip.text = element_text(colour = 'black')) +   theme(strip.background =element_rect(fill="white"))
p <- p1 + scale_color_manual(values=c("{0,k}-mer"="#CC3300","{1,k-1}-mer"="#6BAED6","{2,k-2}-mer"="#4292C6","{3,k-3}-mer"="#2171B5","{4,k-4}-mer"="#08519C","{5,k-5}-mer"="#08306B"))
db2 <- ggpar(p, font.x = c(10,"black"),font.y = c(10,"black"),font.y.tickslab = c(4,"black"), legend = "none")
db2

png("AUCS_dataset2.png", width = 2000, height = 2000, res = 600)
db2
dev.off()

### plots for datasetC
df <- fread(file="AUCS_datasetC.csv",header=T)
df$frag_t = factor(df$frag, levels=c('300 bp','1.000 bp','3.000 bp','5.000 bp','10.000 bp'))

p2 <- ggplot(df, aes(x=frag_t, y=mean, group=type, color=type)) +  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, position=position_dodge(0.05)) + 
  geom_line() + geom_point() + theme_linedraw() + facet_grid(k ~ .) + theme(strip.text.y.left = element_text(angle = 0))
p3 <- p2 + labs(x=NULL, y = "mean AUC", fill = NULL) +  theme(axis.text = element_text(colour = "black", size=11,family='Arial')) 
p5 <- p3 + theme(axis.title = element_text(family = "Arial", size = 11), axis.text = element_text(family = "Arial"), axis.text.x = element_text(family = "Arial"), legend.text = element_text(size = 11, family = "Arial"),  strip.text.y = element_text(size=rel(3.5)))+
  scale_color_manual(values=c("{0,k}-mer"="#1B9E77","{1,k-1}-mer"="#9ECAE1","{2,k-2}-mer"="#4292C6","{3,k-3}-mer"="#084594","{4,k-4}-mer"="#DEEBF7","{5,k-5}-mer"="#C6DBEF"))
p1 <- p3 + theme_linedraw() + theme(strip.text = element_text(colour = 'black')) +   theme(strip.background =element_rect(fill="white"))
p <- p1 + scale_color_manual(values=c("{0,k}-mer"="#CC3300","{1,k-1}-mer"="#6BAED6","{2,k-2}-mer"="#4292C6","{3,k-3}-mer"="#2171B5","{4,k-4}-mer"="#08519C","{5,k-5}-mer"="#08306B"))
db3 <- ggpar(p, font.x = c(10,"black"),font.y = c(10,"black"),font.y.tickslab = c(4,"black"), legend = "none")
db3

png("AUCS_dataset3.png", width = 2000, height = 2000, res = 600)
db3
dev.off()

###

q<- ggarrange(db1, db2, db3, labels = c("A", "B", "C"),ncol = 1, nrow = 3, heights = c(0.62, 1, 1.1))
q

png("Figure2_Andrade_et_al.png", width = 2000, height = 2500, res = 400)
q
dev.off()

