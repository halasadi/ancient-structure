
### SQUAREM implementation (Main function)

#setwd('/Users/kushal/Documents/ancient-structure/src')
library(gtools)
library(SQUAREM)
source('update_EM.R')
source('update_squarem.R')
source('simulate_model.R')
source('loglik_em.R')
source('simplex_functions.R')

npop=5;
nsamp_per_pop=50;
nclusters = 3
maxscale = 10
# omega is defined to be the admixture proportions
# size = nsamp * nclusters

omega <- matrix(data=0,nrow=0,ncol=nclusters);
for(pop in 1:npop) {
omega <- rbind(omega,rdirichlet(nsamp_per_pop,sample(1:maxscale, replace=T, nclusters)));
}

simulate_allele_freq <- function(alpha, nSNPs){
  #mu = 1.25e-8
  #alpha = 4*Ne*mu
  beta = alpha
  return(rbeta(nSNPs, shape1 = alpha, shape2 = beta))
}

# we're working with genotype data so equal number of variants across freq range
alpha <- c(1,1,1,1);

# size = nclusters x nSNPs
freq_mat <- t(matrix(unlist(lapply(1:nclusters, function(n) simulate_allele_freq(alpha[n],nSNPs = 200))),ncol=nclusters));


data <- simulate_binomial_model(omega,freq_mat)


ancient_structure <- function(geno_data, f_known, K_unknown, max_iter)
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
  param_vec_in <- c(as.vector(rev_q_initial),as.vector(f_unknown_initial),as.vector(f_known));
  
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
  f_unknown <- matrix(temp[0:(nSNPs*K_unknown)], nrow=nSNPs, ncol=K_unknown)
  beg = nSNPs*K_unknown + 1
  end = nSNPs*K_pooled
  f_known <- matrix(temp[beg:end], nrow=nSNPs, ncol=K_known)
  
  out <- list("q"=q, "f_known"=f_known,"f_unknown"=f_unknown);
  
  return(out)
  
}

K_known = 1
K_unknown = 2;
K_pooled <- K_known +K_unknown;
out <- ancient_structure(geno_data = data, f_known = t(freq_mat[(K_unknown+1):K_pooled,]), K_unknown = K_unknown, max_iter = 500)

barplot(t(omega),col=2:(nclusters+1),axisnames=F,space=0,border=NA,main=paste("true structure: No. of clusters=",nclusters),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
barplot(t(out$q),col=2:(nclusters+1),axisnames=F,space=0,border=NA,main=paste("estd structure: No. of clusters=",nclusters),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
