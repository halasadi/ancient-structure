
library(data.table)
data.clst <- read.csv('../exploratory/data.clst', stringsAsFactors = FALSE, header=FALSE, sep = "\t");
data.fam <- read.table('../exploratory/haak_fig3_set.fam');

pop_ids <- data.clst[match(data.fam$V2, data.clst$V2),3]

Q_data <- as.matrix(read.table("../internal_data/q_known_3_unknown_2.txt"));
#Q_data <- as.matrix(read.table("../exploratory/Admixture/haak_fig3_set.5.Q"))

K=5
barplot(t(Q_data),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=0.5,cex.main=1.4);

labels = match(unique(pop_ids), pop_ids);
abline(v=labels-1)

labels_low=labels-1;
labels_up=c(labels_low[2:length(labels_low)],dim(Q_data)[1]);
mid_point <- labels_low +0.5*(labels_up-labels_low);
axis(1,at=mid_point, unique(pop_ids),las=2,cex.axis=0.5);

Q_data_pop <- sapply(1:K, function(num) tapply(Q_data[,num], pop_ids,mean));

ancient_pop <- unique(pop_ids)[24:44];
modern_pop <- unique(pop_ids)[1:23];

Q_data_pop_ancient <- Q_data_pop[match(ancient_pop,rownames(Q_data_pop)),];
Q_data_pop_modern <-  Q_data_pop[-match(ancient_pop,rownames(Q_data_pop)),];

barplot(t(Q_data_pop_ancient),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=0.5,cex.main=1.4);
abline(v=1:length(ancient_pop))
mid_point =seq(0.5,dim(Q_data_pop_ancient)[1]-0.5,length.out=dim(Q_data_pop_ancient)[1]);
axis(1,at=mid_point, ancient_pop,las=2,cex.axis=0.5);

barplot(t(Q_data_pop_modern),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=0.5,cex.main=1.4);
abline(v=1:length(modern_pop))
mid_point =seq(0.5,dim(Q_data_pop_modern)[1]-0.5,length.out=dim(Q_data_pop_modern)[1]);
axis(1,at=mid_point, modern_pop,las=2,cex.axis=0.5);

#### compare with the plot for admixture with k=3

P_data <- data.frame(fread(paste0('../exploratory/Admixture/haak_fig3_set.',K,'.P')));
Q_data <- data.frame(fread(paste0('../exploratory/Admixture/haak_fig3_set.',K,'.Q')));

K=3
barplot(t(Q_data),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=0.5,cex.main=1.4);

labels = match(unique(pop_ids), pop_ids);
abline(v=labels-1)

labels_low=labels-1;
labels_up=c(labels_low[2:length(labels_low)],dim(Q_data)[1]);
mid_point <- labels_low +0.5*(labels_up-labels_low);
axis(1,at=mid_point, unique(pop_ids),las=2,cex.axis=0.5);

Q_data_pop <- sapply(1:3, function(num) tapply(Q_data[,num], pop_ids,mean));

ancient_pop <- unique(pop_ids)[24:44];
modern_pop <- unique(pop_ids)[1:23];

Q_data_pop_ancient <- Q_data_pop[match(ancient_pop,rownames(Q_data_pop)),];
Q_data_pop_modern <-  Q_data_pop[-match(ancient_pop,rownames(Q_data_pop)),];

barplot(t(Q_data_pop_ancient),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=0.5,cex.main=1.4);
abline(v=1:length(ancient_pop))
mid_point =seq(0.5,dim(Q_data_pop_ancient)[1]-0.5,length.out=dim(Q_data_pop_ancient)[1]);
axis(1,at=mid_point, ancient_pop,las=2,cex.axis=0.5);

barplot(t(Q_data_pop_modern),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",K),las=1,ylim=c(0,1),cex.axis=0.5,cex.main=1.4);
abline(v=1:length(modern_pop))
mid_point =seq(0.5,dim(Q_data_pop_modern)[1]-0.5,length.out=dim(Q_data_pop_modern)[1]);
axis(1,at=mid_point, modern_pop,las=2,cex.axis=0.5);


