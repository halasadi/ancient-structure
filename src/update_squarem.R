
library(SQUAREM)

# as vector goes column wise
# param_vec_in <- c(as.vector(q_in),as.vector(f_unknown_in),as.vector(f_known_in));

update_squarem <- function(param_vec_in, geno_data, nsamp = 250, K_unknown = 0, K_known = 4, nSNPs = 5)
{
  K_pooled <- K_known +K_unknown;
  q_in = matrix(param_vec_in[(1:(nsamp*nSNPs))],nrow = nsamp, ncol = (K_pooled));
  temp <- param_vec_in[-(1:(nsamp*nSNPs))];
  f_unknown_in <- matrix(temp[1:(nSNPs*K_unknown)], nrow=nSNPs, ncol=K_unknown)
  f_known_in <- matrix(temp[-(1:(nSNPs*K_unknown))], nrow=nSNPs, ncol=K_known)
  out <- update_EM(q_in, f_unknown_in,f_known_in,geno_data);
  param_vec_out <- c(as.vector(out$q),as.vector(out$f_unknown),as.vector(out$f_known));
  return(param_vec_out)
}

