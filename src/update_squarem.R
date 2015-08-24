
library(SQUAREM)

# as vector goes column wise
# param_vec_in <- c(as.vector(q_in),as.vector(f_unknown_in),as.vector(f_known_in));

# SQUAREM wrapper function

# @param_vec_in : a vector of size nsamp*n_clus+n_clus_known*n_genes+n_clus_unknown*ngenes: rev_trans(q), logit(f_known), logit(f_unknown) vectorized
# update_EM is the main function 
# This function prepares the input for update_EM and then reverse transforms the output for squarem interpolation



update_squarem <- function(param_vec_in, geno_data, nsamp, K_unknown, K_known, nSNPs)
{
  
  K_pooled <- K_known +K_unknown;
  
  # rev_q_in is nsamp*(K_pooled-1) matrix, each row is a reverse transform on q (topic prop matrix)
  # so has one less co-ordinate as it is not constrained
  
  rev_q_in = matrix(param_vec_in[(1:(nsamp*(K_pooled-1)))],nrow = nsamp, ncol = (K_pooled-1));
  
  # transforming the rev_q_in to q_in by adding the sum=1 constraint
  
  q_in <- t(apply(rev_q_in, 1,function(x) transform(x)));
  
  # extracting the unknown allele frequencies data from the param_vec_in object and inv transforming it
  
  temp <- param_vec_in[-(1:(nsamp*(K_pooled-1)))];
  f_unknown_in <- inv.logit(matrix(temp[0:(nSNPs*K_unknown)], nrow=nSNPs, ncol=K_unknown));
  
  # extracting the known allele frequencies data from the param_vec_in object and inv transforming it
  
  beg = nSNPs*K_unknown + 1;
  end = nSNPs*K_pooled;
  f_known_in <- inv.logit(matrix(temp[beg:end], nrow=nSNPs, ncol=K_known));
  
  # using the EM update scheme- the main function 
  out <- update_EM(q_in, f_unknown_in,f_known_in,geno_data);
  rev_q_out <- as.matrix(t(apply(out$q, 1, function(x) reverse_transform(x))));
  param_vec_out <- c(as.vector(rev_q_out),as.vector(logit(out$f_unknown)),as.vector(logit(out$f_known)));
  return(param_vec_out)
}

