
#### Heterozygosity plot across chromosomes

library(data.table)
all.hwe <- read.table('all.hwe',header=TRUE)

snp_locations <- unlist(lapply(all.hwe$SNP, function(x) as.numeric(strsplit(as.character(x),":")[[1]][2])))

data <- cbind.data.frame(all.hwe$CHR,snp_locations,all.hwe$P);
colnames(data)=c("CHR","BP","P");

library(qqman)

manhattan(data)

data2 <- cbind.data.frame(all.hwe$CHR,snp_locations,all.hwe$O.HET.);
colnames(data2)=c("CHR","BP","P");

manhattan(data2, suggestiveline = FALSE,genomewideline = FALSE,logp = FALSE,
          ylab="Obs. het", ylim=c(0,0.6))


