library(gtools)
library(SQUAREM)
library(parallel)
library(boot)
library(plyr)
source('update_EM.R')
source('update_squarem.R')
source('simulate_model.R')
source('loglik_em.R')
source('simplex_functions.R')

ancient_structure <- function(geno_data, f_known, K_unknown, max_iter, eps=1e-04)
{
  f_known <- matrix(f_known, nrow=dim(geno_data)[2])
  K_known = dim(f_known)[2]
  K_pooled = K_known + K_unknown
  nsamp <- dim(geno_data)[1];
  nSNPs <- dim(geno_data)[2];
  f_unknown_initial = matrix(nrow = nSNPs, ncol = K_unknown, runif(nSNPs*K_unknown))
  q_initial = matrix(nrow = nsamp, ncol = K_pooled, rdirichlet(nsamp, rep(1/K_pooled,K_pooled)));
  rev_q_initial <-as.matrix(t(apply(q_initial, 1, function(x) reverse_transform(x))));
  # column wise
  ind0 <- which(f_known < eps, arr.ind=T);
  ind1 <- which(f_known > (1-eps), arr.ind=T);
  f_known[ind1] = 1-eps
  f_known[ind0] = eps
  logit_f_known <- logit(f_known);
  logit_f_unknown_initial <- logit(f_unknown_initial);
  param_vec_in <- c(as.vector(rev_q_initial),as.vector(logit_f_unknown_initial),as.vector(logit_f_known));
  
  
  res <- squarem(par=as.numeric(param_vec_in),
                 fixptfn=update_squarem,
                 #objfn= loglik_squarem,
                 geno_data = geno_data,
                 nsamp = nsamp,
                 nSNPs = nSNPs,
                 K_unknown = K_unknown,
                 K_known = K_known,
                 control=list(maxiter = max_iter, trace = FALSE));
  
  ## TESTING ##
  
  #  for(iter in 1:max_iter) { 
  #    out = update_squarem(param_vec_in, geno_data, nsamp, K_unknown, K_known, nSNPs);
  #    param_vec_in <- out; 
  #    print(iter);
  #    rev_q_temp = matrix(out[(1:(nsamp*(K_pooled-1)))],nrow = nsamp, ncol = (K_pooled-1));
  #    q_temp <- t(apply(rev_q_temp, 1,function(x) transform(x)));
  #    print(q_temp)}
  
  ## END TESTING ##
  
  rev_q = matrix(res$par[(1:(nsamp*(K_pooled-1)))],nrow = nsamp, ncol = (K_pooled-1));
  q <- t(apply(rev_q, 1,function(x) transform(x)));
  # q = matrix(res$par[(1:(nsamp*K_pooled))],nrow = nsamp, ncol = K_pooled);
  temp <- res$par[-(1:(nsamp*(K_pooled-1)))];
  f_unknown <- inv.logit(matrix(temp[0:(nSNPs*K_unknown)], nrow=nSNPs, ncol=K_unknown))
  beg = nSNPs*K_unknown + 1;
  end = nSNPs*K_pooled;
  f_known <- inv.logit(matrix(temp[beg:end], nrow=nSNPs, ncol=K_known))
  
  out <- list("q"=q, "f_known"=f_known,"f_unknown"=f_unknown);
  
  return(out)
  
}
