
library(SQUAREM)

# as vector goes column wise
# param_vec_in <- c(as.vector(q_in),as.vector(f_unknown_in),as.vector(f_known_in));

update_squarem <- function(param_vec_in, geno_data, nsamp = 250, K_unknown = 0, K_known = 4, nSNPs = 5)
{
  K_pooled <- K_known +K_unknown;
  rev_q_in = matrix(param_vec_in[(1:(nsamp*(K_pooled-1)))],nrow = nsamp, ncol = (K_pooled-1));
  q_in <- t(apply(rev_q_in, 1,function(x) transform(x)));
  temp <- param_vec_in[-(1:(nsamp*(K_pooled-1)))];
  f_unknown_in <- inv.logit(matrix(temp[0:(nSNPs*K_unknown)], nrow=nSNPs, ncol=K_unknown));
  beg = nSNPs*K_unknown + 1;
  end = nSNPs*K_pooled;
  f_known_in <- inv.logit(matrix(temp[beg:end], nrow=nSNPs, ncol=K_known));
  out <- update_EM(q_in, f_unknown_in,f_known_in,geno_data);
  rev_q_out <- as.matrix(t(apply(out$q, 1, function(x) reverse_transform(x))));
  param_vec_out <- c(as.vector(rev_q_out),as.vector(logit(out$f_unknown)),as.vector(logit(out$f_known)));
  return(param_vec_out)
}

