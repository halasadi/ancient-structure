

library(gtools)
library(SQUAREM)
library(parallel)
library(boot)
library(plyr)
set.seed(10)
source('utilities.R')
# The main function that fits the ancient structure and gives the final topic proportions and allele frequencies

ancient_structure <- function(geno_data, K_unknown, pop_labs, source_labs, max_iter, eps=1e-04, use_squarem=FALSE)
{
  K_known <- length(source_labs);
  K_pooled = K_known + K_unknown
  nsamp <- dim(geno_data)[1];
  nSNPs <- dim(geno_data)[2];

  # determining initial values of f_unknown and q
  f_pooled_initial = matrix(nrow = nSNPs, ncol = K_pooled, runif(nSNPs*K_unknown))
  q_initial = matrix(nrow = nsamp, ncol = K_pooled, rdirichlet(nsamp, rep(1/K_pooled,K_pooled)));
  q_initial <- fix_clus_mem(q_initial,pop_labs, source_labs);
  rev_q_initial <-matrix(t(apply(q_initial, 1, function(x) reverse_transform(x))), nrow=nsamp);
  
  # transforming f_known and f_unknown to logit form to make it unconstrained, needed for squarem input
  logit_f_pooled_initial <- logit(f_pooled_initial);

  # pooling the transformed q and the logit transformed f_known and f_unknown

  param_vec_in <- c(as.vector(rev_q_initial),as.vector(logit_f_pooled_initial));

  # using squarem

  if(use_squarem)
  {
    res <- squarem(par=as.numeric(param_vec_in),
                   fixptfn=update_squarem,
                   objfn= loglik_squarem,
                   geno_data = geno_data,
                   nsamp = nsamp,
                   nSNPs = nSNPs,
                   K_unknown = K_unknown,
                   pop_labs =pop_labs,
                   source_labs=source_labs,
                   control=list(maxiter = max_iter, trace = FALSE, square=FALSE, tol=1e-10));
    rev_q = matrix(res$par[(1:(nsamp*(K_pooled-1)))],nrow = nsamp, ncol = (K_pooled-1));
    q <- t(apply(rev_q, 1,function(x) transform(x)));
    # q = matrix(res$par[(1:(nsamp*K_pooled))],nrow = nsamp, ncol = K_pooled);
    temp <- res$par[-(1:(nsamp*(K_pooled-1)))];
    f_pooled <- inv.logit(matrix(temp[0:(nSNPs*K_pooled)], nrow=nSNPs, ncol=K_pooled))
    
    out <- list("q"=q, "f_pooled"=f_pooled);

    return(out)
  }

  if(!use_squarem)
  {
    for(iter in 1:max_iter) {
      out = update_squarem(param_vec_in, geno_data, pop_labs, source_labs, K_unknown);
      loglik <- loglik_squarem(param_vec_in, geno_data, pop_labs, source_labs, K_unknown);
      cat('Neg Loglikelihood at iteration',iter,':',loglik,'\n');
      param_vec_in <- out;
      }

      rev_q = matrix(param_vec_in[(1:(nsamp*(K_pooled-1)))],nrow = nsamp, ncol = (K_pooled-1));
      q <- t(apply(rev_q, 1,function(x) transform(x)));
      temp <- param_vec_in[-(1:(nsamp*(K_pooled-1)))];
      f_pooled <- inv.logit(matrix(temp[0:(nSNPs*K_pooled)], nrow=nSNPs, ncol=K_pooled))
      outlist <- list("q"=q, "f_pooled"=f_pooled);
      return(outlist)
  }



}
