
#loglik_em <- function(q_in, f_unknown_in, f_known_in, geno_data)
#{
#  f_pooled_in <- cbind(f_unknown_in, f_known_in);
#  prod <- q_in %*% f_pooled_in;
#  return( sum(log(prod)*geno_data + (2-geno_data)*log(1-prod)) )
#}

# @param_vec_in : a vector of size nsamp*n_clus+n_clus_known*n_genes+n_clus_unknown*ngenes: rev_trans(q), logit(f_known), logit(f_unknown) vectorized

loglik_squarem <- function(param_vec_in, geno_data, nsamp, K_unknown, K_known, nSNPs)
{
  K_pooled <- K_known + K_unknown;
  
  # rev_q_in is nsamp*(K_pooled-1) matrix, each row is a reverse transform on q (topic prop matrix)
  # so has one less co-ordinate as it is not constrained
  
  rev_q_in = matrix(param_vec_in[(1:(nsamp*(K_pooled-1)))],nrow = nsamp, ncol = (K_pooled-1));
  
 
  # transforming the rev_q_in to q_in by adding the sum=1 constraint
  
  q_in <- t(apply(rev_q_in, 1,function(x) transform(x)));
  
  
  # extracting the unknown allele frequencies data from the param_vec_in object
  
  temp <- param_vec_in[-(1:(nsamp*(K_pooled-1)))];
  f_unknown_in <- matrix(temp[0:(nSNPs*K_unknown)], nrow=nSNPs, ncol=K_unknown)
  
  
  # extracting the known allele frequencies data from the param_vec_in object
  
  beg = nSNPs*K_unknown + 1
  end = nSNPs*K_pooled
  f_known_in <- matrix(temp[beg:end], nrow=nSNPs, ncol=K_known)
  
  # pooling the unknown and known allele frequencies
  
  f_pooled_in <- cbind(f_unknown_in, f_known_in);
  
  # loglikelihood computation
  prod <- q_in %*% t(f_pooled_in);
  return( sum(log(prod)*geno_data + (2-geno_data)*log(1-prod)))
}