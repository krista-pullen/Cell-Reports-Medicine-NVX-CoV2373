library(viridis)
library(gplots)
library(robustbase)
library(dplyr)
library(impute)
library(ggplot2)
library(reshape2)

### PREPROCESS DATA ###

full <- read.csv(file ='/Users/kpullen/novovax-human-data.csv',header= TRUE)
full<-full[c(1:237),]
peak<-full[full$TP=='2',]

#remove variant data
wt<-select(peak,-c(contains("E406Q"),contains("N440K"),contains("d69"),contains("E484K"),contains("HKU"),contains("EBOV"),contains("SA"),contains("N501Y"),contains("D614G"),contains("F486A"),contains("K417"),contains("F490K"),contains("Q493R"),contains("N487R"),contains("NVX")))
data<-wt[,c(8:43)]

#log transform and normalize data
data[,c(1:30,33)]<-log10(data[,c(1:30,33)])
data<-data[,c(31:36,1:30)]
X<-scale(data, center = TRUE, scale = TRUE)
X<-data.frame(X)

# Reorder and label features with common names
X<-X[,c(1:6,9:10,7,11,8,34:35,32,35,33,29:30,27,31,28,14:15,12,16,13,19:20,17,21,18,24:25,22,26,23)]

type<-rep(c("Functions","IgG1","IgG2","IgG3","FcgR2a","FcgR2b","FcgR3a"), c(6,5,5,5,5,5,5))

rc <-rep(c('#AF93C7','#704793',"#35085B"), c(length(wt$Arm[wt$Arm=='ARM B']),length(wt$Arm[wt$Arm=='ARM D']),length(wt$Arm[wt$Arm=='ARM C'])))

### FIGURE 6B ###

heatmap.2(data.matrix(X),dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',keysize=1,
          col = magma(n=15,begin = 0, end = 1),margins = c(5,28),RowSideColors=rc,breaks = seq(-3, 7,length.out = 16))

### FIGURE 6D ###

X_percent<-X
for (ind_feat in 1:ncol(X_percent)) {
  X_percent[, ind_feat] <- percent_rank(X_percent[, ind_feat])
}

group25<-X_percent[wt$Arm=='ARM B',]
group25M<-X_percent[wt$Arm=='ARM D',]
group5M<-X_percent[wt$Arm=='ARM C',]

group25<-melt(group25)
group25$variable <- factor(group25$variable, levels = unique(group25$variable))
group25$value <- as.numeric(group25$value)
group25<- aggregate(group25$value, by=list(group25$variable), FUN=mean)
colnames(group25)<-c('variable','value')
ggplot(group25, aes(factor(variable), value, fill = factor(variable))) +geom_bar(stat='identity',width = 1)+coord_polar() + geom_hline(yintercept = seq(0,0.8, by = 0.2),color = 'grey', size = 1) +geom_vline(xintercept = seq(.5, 20, by = 1),color = 'grey', size = 1) +theme(panel.background = element_blank())

group25M<-melt(group25M)
group25M$variable <- factor(group25M$variable, levels = unique(group25M$variable))
group25M$value <- as.numeric(group25M$value)
group25M<- aggregate(group25M$value, by=list(group25M$variable), FUN=mean)
colnames(group25M)<-c('variable','value')


ggplot(group25M, aes(factor(variable), value, fill = factor(variable))) +geom_bar(stat='identity',width = 1)+coord_polar() + geom_hline(yintercept = seq(0,0.8, by = 0.2),color = 'grey', size = 1) +geom_vline(xintercept = seq(.5, 20, by = 1),color = 'grey', size = 1) +theme(panel.background = element_blank())

group5M<-melt(group5M)
group5M$variable <- factor(group5M$variable, levels = unique(group5M$variable))
group5M$value <- as.numeric(group5M$value)
group5M<- aggregate(group5M$value, by=list(group5M$variable), FUN=mean)
colnames(group5M)<-c('variable','value')

group5M<-group5M[-c(1:6),]
ggplot(group5M, aes(factor(variable), value, fill = factor(variable))) +geom_bar(stat='identity',width = 1)+coord_polar() + geom_hline(yintercept = seq(0,0.8, by = 0.2),color = 'grey', size = 1) +geom_vline(xintercept = seq(.5, 20, by = 1),color = 'grey', size = 1) +theme(panel.background = element_blank())


time<-full[full$Arm == 'ARM C',]

wt<-select(time,-c(contains("E406Q"),contains("N440K"),contains("d69"),contains("E484K"),contains("HKU"),contains("EBOV"),contains("SA"),contains("N501Y"),contains("D614G"),contains("F486A"),contains("K417"),contains("F490K"),contains("Q493R"),contains("N487R"),contains("NVX")))
wt<-wt[,c(8:43)]
wt[,c(1:30,33)]<-log10(wt[,c(1:30,33)])
wt<-wt[,c(31:36,1:30)]
X<-scale(wt, center = TRUE, scale = TRUE)
X<-data.frame(X)
X<-X[,c(1:6,9:10,7,11,8,34:35,32,35,33,29:30,27,31,28,14:15,12,16,13,19:20,17,21,18,24:25,22,26,23)]

X_percent<-X
for (ind_feat in 1:ncol(X_percent)) {
  X_percent[, ind_feat] <- percent_rank(X_percent[, ind_feat])
}

group1<-X_percent[time$TP=='1',]
group2<-X_percent[time$TP=='2',]
group3<-X_percent[time$TP=='3',]

group1<-melt(group1)
group1$variable <- factor(group1$variable, levels = unique(group1$variable))
group1$value <- as.numeric(group1$value)
group1<- aggregate(group1$value, by=list(group1$variable), FUN=mean)
colnames(group1)<-c('variable','value')

ggplot(group1, aes(factor(variable), value, fill = factor(variable))) +geom_bar(stat='identity',width = 1)+coord_polar() + geom_hline(yintercept = seq(0,0.8, by = 0.2),color = 'grey', size = 1) +geom_vline(xintercept = seq(.5, 20, by = 1),color = 'grey', size = 1) +theme(panel.background = element_blank())

group2<-melt(group2)
group2$variable <- factor(group2$variable, levels = unique(group2$variable))
group2$value <- as.numeric(group2$value)
group2<- aggregate(group2$value, by=list(group2$variable), FUN=mean)
colnames(group2)<-c('variable','value')

ggplot(group2, aes(factor(variable), value, fill = factor(variable))) +geom_bar(stat='identity',width = 1)+coord_polar() + geom_hline(yintercept = seq(0,0.8, by = 0.2),color = 'grey', size = 1) +geom_vline(xintercept = seq(.5, 20, by = 1),color = 'grey', size = 1) +theme(panel.background = element_blank())

group3<-melt(group3)
group3$variable <- factor(group3$variable, levels = unique(group3$variable))
group3$value <- as.numeric(group3$value)
group3<- aggregate(group3$value, by=list(group3$variable), FUN=mean)
colnames(group3)<-c('variable','value')
#group5M<-group5M[-c(1:6),]
ggplot(group3, aes(factor(variable), value, fill = factor(variable))) +geom_bar(stat='identity',width = 1)+coord_polar() + geom_hline(yintercept = seq(0,0.8, by = 0.2),color = 'grey', size = 1) +geom_vline(xintercept = seq(.5, 20, by = 1),color = 'grey', size = 1) +theme(panel.background = element_blank())

### FIGURE 6D ###

time<-full[full$Arm == 'ARM C',]

wt<-select(time,-c(contains("E406Q"),contains("N440K"),contains("d69"),contains("E484K"),contains("HKU"),contains("EBOV"),contains("SA"),contains("N501Y"),contains("D614G"),contains("F486A"),contains("K417"),contains("F490K"),contains("Q493R"),contains("N487R"),contains("NVX")))
wt<-wt[,c(8:43)]
wt[,c(1:30,33)]<-log10(wt[,c(1:30,33)])
wt<-wt[,c(31:36,1:30)]
X<-scale(wt, center = TRUE, scale = TRUE)
X<-data.frame(X)
X<-X[,c(1:6,9:10,7,11,8,34:35,32,35,33,29:30,27,31,28,14:15,12,16,13,19:20,17,21,18,24:25,22,26,23)]

X_percent<-X
for (ind_feat in 1:ncol(X_percent)) {
  X_percent[, ind_feat] <- percent_rank(X_percent[, ind_feat])
}

group1<-X_percent[time$TP=='1',]
group2<-X_percent[time$TP=='2',]
group3<-X_percent[time$TP=='3',]

group1<-melt(group1)
group1$variable <- factor(group1$variable, levels = unique(group1$variable))
group1$value <- as.numeric(group1$value)
group1<- aggregate(group1$value, by=list(group1$variable), FUN=mean)
colnames(group1)<-c('variable','value')


ggplot(group1, aes(factor(variable), value, fill = factor(variable))) +geom_bar(stat='identity',width = 1)+coord_polar() + geom_hline(yintercept = seq(0,0.8, by = 0.2),color = 'grey', size = 1) +geom_vline(xintercept = seq(.5, 20, by = 1),color = 'grey', size = 1) +theme(panel.background = element_blank())

group2<-melt(group2)
group2$variable <- factor(group2$variable, levels = unique(group2$variable))
group2$value <- as.numeric(group2$value)
group2<- aggregate(group2$value, by=list(group2$variable), FUN=mean)
colnames(group2)<-c('variable','value')


ggplot(group2, aes(factor(variable), value, fill = factor(variable))) +geom_bar(stat='identity',width = 1)+coord_polar() + geom_hline(yintercept = seq(0,0.8, by = 0.2),color = 'grey', size = 1) +geom_vline(xintercept = seq(.5, 20, by = 1),color = 'grey', size = 1) +theme(panel.background = element_blank())

group3<-melt(group3)
group3$variable <- factor(group3$variable, levels = unique(group3$variable))
group3$value <- as.numeric(group3$value)
group3<- aggregate(group3$value, by=list(group3$variable), FUN=mean)
colnames(group3)<-c('variable','value')

ggplot(group3, aes(factor(variable), value, fill = factor(variable))) +geom_bar(stat='identity',width = 1)+coord_polar() + geom_hline(yintercept = seq(0,0.8, by = 0.2),color = 'grey', size = 1) +geom_vline(xintercept = seq(.5, 20, by = 1),color = 'grey', size = 1) +theme(panel.background = element_blank())


