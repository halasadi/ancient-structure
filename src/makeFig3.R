library(gtools)
library(SQUAREM)
library(parallel)
library(data.table)
#source('update_EM.R')
#source('update_squarem.R')
#source('loglik_em.R')
#source('simplex_functions.R')
source('ancient_structure.R')

data <- as.matrix(data.frame(fread("../external_data/imputed_haak_fig3_genotypes.txt"),row.names=TRUE));
#data <- data[, 1:10000];
freq_mat <- as.matrix(data.frame(fread("../external_data/ancestral_freqs.txt"), row.names = TRUE));
#freq_mat <- freq_mat[,1:10000]

removeDup <- function(M){
    inter0 = (M[1,] == 0)
    inter = (M[1,] == 1)
    for (i in 2:nrow(M)){
        inter0 = inter0 * (M[i,] == 0)
        inter = inter & (M[i,] == 1)
    }
    return(as.numeric(unique(c(which(inter == TRUE), which(inter0 == TRUE)))))
}

snpsToRemove = removeDup(freq_mat)

data <- data[,-snpsToRemove]
freq_mat <- freq_mat[,-snpsToRemove]
data <- data[,1:20000]
freq_mat <- freq_mat[1:2,1:20000]

K_unknown = 0;

out <- ancient_structure(geno_data = data, f_known = t(freq_mat), K_unknown = K_unknown, max_iter = 100, eps=1e-02, use_squarem=FALSE)

#write.table(out$q, file="q.txt");
write.table(out$q, file = "../internal_data/q_one_unknown.txt")
write.table(out$q, file = "../internal_data/q_three_unknown.txt")

#nclusters <- 3

#barplot(t(out$q),col=2:(nclusters+1),axisnames=F,space=0,border=NA,main=paste("estd structure: No. of clusters=",nclusters),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


