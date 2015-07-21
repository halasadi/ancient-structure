
### SQUAREM implementation (Main function)

library(gtools)
library(SQUAREM)

ancient_structure <- function(geno_data, f_known, K_unknown, max_iter)
{
  K_known = dim(f_known)[2]
  K_pooled = K_known + K_unknown
  nsamp <- dim(geno_data)[1];
  nSNPs <- dim(geno_data)[2];
  f_unknown_initial = matrix(nrow = nSNPs, ncol = K_unknown, runif(nSNPs*K_unknown))
  q_initial = matrix(nrow = nsamp, ncol = K_pooled, rdirichlet(nsamp, rep(1/K_pooled,K_pooled)));
  param_vec_in <- c(as.vector(q_initial),as.vector(f_unknown_initial),as.vector(f_known));
  
  res <- squarem(par=as.numeric(param_vec_in),
                             fixptfn=update_squarem,
                             objfn= loglik_squarem,
                             geno_data = geno_data,
                             nsamp = nsamp,
                             nSNPs = nSNPs,
                             K_unknown = K_unknown,
                             K_known = K_known,
                             control=list(maxiter = max_iter, trace = FALSE));
  q = matrix(res$par[(1:(nsamp*nSNPs))],nrow = nsamp, ncol = nSNPs);
  temp <- res$par[-(1:(nsamp*nSNPs))];
  f_unknown <- matrix(temp[1:(nSNPs*K_unknown)], nrow=nSNPs, ncol=K_unknown)
  f_known <- matrix(temp[-(1:(nSNPs*K_unknown))], nrow=nSNPs, ncol=K_known)
  
  out <- list("q"=q, "f_known"=f_known,"f_unknown"=f_unknown);
  
  return(out)
  
}

ancient_structure(geno_data = geno_data, f_known = f_known, K_unknown = 0, max_iter = 10)
