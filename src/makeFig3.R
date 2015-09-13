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

n_obs <- c(18,24,2);

out <- ancient_structure(geno_data = data, f_obs = t(freq_mat), n_obs = n_obs, K_unknown = K_unknown, max_iter = 2, eps=1e-02, use_squarem=FALSE)

write.table(out$q, file = paste0("../internal_data/sampling_correction/q_known_", K_known,
                       "_unknown_", K_unknown, ".txt"))



