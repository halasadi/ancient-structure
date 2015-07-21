
loglik_em <- function(q_in, f_unknown_in, f_known_in, geno_data)
{
  f_pooled_in <- cbind(f_unknown_in, f_known_in);
  prod <- q_in %*% f_pooled_in;
  return( sum(prod*geno_data +(2-geno_data)*(1-prod)) )
}


loglik_squarem <- function(param_vec_in, geno_data, nsamp = 250, K_unknown = 0, K_known = 4, nSNPs = 5)
{
  K_pooled <- K_known +K_unknown;
  q_in = matrix(param_vec_in[(1:(nsamp*nSNPs))],nrow = nsamp, ncol = K_pooled);
  temp <- param_vec_in[-(1:(nsamp*nSNPs))];
  f_unknown_in <- matrix(temp[1:(nSNPs*K_unknown)], nrow=nSNPs, ncol=K_unknown)
  f_known_in <- matrix(temp[-(1:(nSNPs*K_unknown))], nrow=nSNPs, ncol=K_known)
  f_pooled_in <- cbind(f_unknown_in, f_known_in);
  prod <- q_in %*% t(f_pooled_in);
  return( sum(prod*geno_data +(2-geno_data)*(1-prod)) )
}