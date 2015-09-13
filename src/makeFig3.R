library(gtools)
library(SQUAREM)
library(parallel)
library(data.table)
source('ancient_structure.R')

args <- commandArgs(trailingOnly = TRUE)
print(args)
K_known = as.integer(args[1])
K_unknown = as.integer(args[2])



data <- as.matrix(data.frame(fread("../external_data/imputed_haak_fig3_genotypes.txt"),row.names=TRUE));

freq_mat <- as.matrix(data.frame(fread("../external_data/ancestral_freqs.txt"), row.names = TRUE));

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
#data <- data[,1:20000]

#freq_mat <- freq_mat[1:K_known,1:10000]
#data <- data[,1:10000];

out <- ancient_structure(geno_data = data, f_known = t(freq_mat), K_unknown = K_unknown, max_iter = 100, eps=1e-02, use_squarem=FALSE)

write.table(out$q, file = paste0("../internal_data/q_known_", K_known,
                       "_unknown_", K_unknown, ".txt"))



