
library(gtools)
library(SQUAREM)
library(parallel)
library(boot)
library(plyr)
source('utilities.R')

# The main function that fits the ancient structure and gives the final topic proportions and allele frequencies

ancient_structure <- function(geno_data, f_obs, n_obs, K_unknown, max_iter, eps=1e-04, use_squarem=FALSE)
{
  f_obs <- matrix(f_obs, nrow=dim(geno_data)[2])
  K_known = dim(f_obs)[2]
  K_pooled = K_known + K_unknown
  nsamp <- dim(geno_data)[1];
  nSNPs <- dim(geno_data)[2];

  # determining initial values of f_unknown and q
  f_unknown_initial = matrix(nrow = nSNPs, ncol = K_unknown, runif(nSNPs*K_unknown));
  q_initial = matrix(nrow = nsamp, ncol = K_pooled, rdirichlet(nsamp, rep(1/K_pooled,K_pooled)));
  rev_q_initial <-as.matrix(t(apply(q_initial, 1, function(x) reverse_transform(x))));

  # trying to add some perturbation 'eps' to counter values of 0 and 1 in allele frequencies (eps user defined)

  ind0 <- which(f_obs < eps, arr.ind=T);
  ind1 <- which(f_obs > (1-eps), arr.ind=T);
  f_obs[ind1] = 1-eps
  f_obs[ind0] = eps

  f_known_initial = f_obs;

  # transforming f_known and f_unknown to logit form to make it unconstrained, needed for squarem input

  logit_f_known_initial <- logit(f_known_initial);
  logit_f_unknown_initial <- logit(f_unknown_initial);

  # pooling the transformed q and the logit transformed f_known and f_unknown

  param_vec_in <- c(as.vector(rev_q_initial),as.vector(logit_f_unknown_initial),as.vector(logit_f_known_initial));

  # using squarem

  if(use_squarem)
  {
    res <- squarem(par=as.numeric(param_vec_in),
                   fixptfn=update_squarem,
                   objfn= loglik_squarem,
                   geno_data = geno_data,
                   f_obs =f_obs,
                   n_obs = n_obs,
                   nsamp = nsamp,
                   nSNPs = nSNPs,
                   K_unknown = K_unknown,
                   K_known = K_known,
                   control=list(maxiter = max_iter, trace = FALSE, square=FALSE, tol=1e-10));
    rev_q = matrix(res$par[(1:(nsamp*(K_pooled-1)))],nrow = nsamp, ncol = (K_pooled-1));
    q <- t(apply(rev_q, 1,function(x) transform(x)));
    # q = matrix(res$par[(1:(nsamp*K_pooled))],nrow = nsamp, ncol = K_pooled);
    temp <- res$par[-(1:(nsamp*(K_pooled-1)))];
    f_unknown <- inv.logit(matrix(temp[0:(nSNPs*K_unknown)], nrow=nSNPs, ncol=K_unknown))
    beg = nSNPs*K_unknown + 1;
    end = nSNPs*K_pooled;
    if(beg >end)
      f_known <- inv.logit(matrix(0, nrow=nSNPs, ncol=K_known));
    if(end >beg)
      f_known <- inv.logit(matrix(temp[beg:end], nrow=nSNPs, ncol=K_known));


    out <- list("q"=q, "f_known"=f_known,"f_unknown"=f_unknown);

    return(out)
  }

  if(!use_squarem)
  {
    for(iter in 1:max_iter) {
      out = update_squarem(param_vec_in, geno_data, f_obs, n_obs, nsamp, K_unknown, K_known, nSNPs);
      loglik <- loglik_squarem(param_vec_in, geno_data, f_obs, n_obs, nsamp, K_unknown, K_known, nSNPs);
      cat('Neg Loglikelihood at iteration',iter,':',loglik,'\n');
      param_vec_in <- out;
      }

      rev_q = matrix(param_vec_in[(1:(nsamp*(K_pooled-1)))],nrow = nsamp, ncol = (K_pooled-1));
      q <- t(apply(rev_q, 1,function(x) transform(x)));
      temp <- param_vec_in[-(1:(nsamp*(K_pooled-1)))];
      f_unknown <- inv.logit(matrix(temp[0:(nSNPs*K_unknown)], nrow=nSNPs, ncol=K_unknown))
      beg = nSNPs*K_unknown + 1;
      end = nSNPs*K_pooled;
      if(beg >end)
        f_known <- inv.logit(matrix(0, nrow=nSNPs, ncol=K_known));
      if(end >beg)
        f_known <- inv.logit(matrix(temp[beg:end], nrow=nSNPs, ncol=K_known));

      out_list <- list("q"=q, "f_known"=f_known,"f_unknown"=f_unknown);

      return(out_list)
  }



}
