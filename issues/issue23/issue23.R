
## Issue 23
setwd('~/Documents/ancient-structure/issues/issue23/')
library(data.table)
data <- as.matrix(data.frame(fread("../../external_data/imputed_haak_fig3_genotypes.txt"),row.names=TRUE));

allele_freq_mat <- as.matrix(data.frame(fread("../../internal_data/pooled_admixture/haak_fig3.LDprune.3.P")));
omega_mat <- as.matrix(data.frame(fread("../../internal_data/pooled_admixture/haak_fig3.LDprune.3.Q")));

geno_data <- data[1:331,];
omega_mat <- omega_mat[1:331,];

loglik_calc <- function(geno_data, omega_mat, allele_freq_mat)
{
  prod <- omega_mat %*% t(allele_freq_mat);
  return(sum(log(prod)*geno_data + (2-geno_data)*log(1-prod)))
}

loglik1 <- loglik_calc(geno_data, omega_mat, allele_freq_mat)
