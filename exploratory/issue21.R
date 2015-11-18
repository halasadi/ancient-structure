
## Build the ancetsral freqs from general admixture plots

setwd("/Users/kushal/Documents/ancient-structure/exploratory/")
library(data.table)
theta <- as.matrix(data.frame(fread("Admixture/haak_fig3.LDprune.3.P")))
omega <- as.matrix(data.frame(fread("Admixture/haak_fig3.LDprune.3.Q")))

data.clst <- read.csv('../external_data/data.clst', stringsAsFactors = FALSE, header=FALSE, sep = "\t");
data.fam <- read.table('../external_data/haak_fig3.LDprune.fam');

pop_ids <- data.clst[match(data.fam$V2, data.clst$V2),3]

## Configuring the ancetsral allele freq of Yamnaya

Yamnaya_freqs <- omega[which(pop_ids=="Yamnaya"),] %*% t(theta)
Yamnaya_freq <- colMeans(Yamnaya_freqs);

## Configuring the ancestral allele freqs of LBK_EN

LBK_freqs <- omega[which(pop_ids=="LBK_EN"),] %*% t(theta)
LBK_freq <- colMeans(LBK_freqs);

## Configuring the ancestral allele freqs of Loschbour

Loschbour_freq <- omega[which(pop_ids=="Loschbour"),] %*% t(theta)


ancestral_freqs <- rbind(Yamnaya_freq, LBK_freq, Loschbour_freq);

write.table(ancestral_freqs, "../exploratory/ancestral_freqs_Yam_LBK_Losch.txt")

## Configuring  the aancestral allele freq of second cluster

green_freq <- theta[,2];

ancestral_freqs <- rbind(Yamnaya_freq, LBK_freq, green_freq)

write.table(ancestral_freqs, "../exploratory/ancestral_freqs_Yam_LBK_green.txt")




