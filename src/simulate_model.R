
# @: function to simulate a N \times L matrix for L sites and N samples 
# day populations and N_a ancient populations. 

npop=5;
nsamp_per_pop=50;
nclusters = 4
maxscale = 10
# omega is defined to be the admixture proportions
# size = nsamp * nclusters
omega = matrix(rbind(rdirichlet(nsamp_per_pop,sample(1:maxscale, replace=T, nclusters)),
                        rdirichlet(nsamp_per_pop,sample(1:maxscale, replace=T, nclusters)),
                        rdirichlet(nsamp_per_pop,sample(1:maxscale, replace=T, nclusters)),
                        rdirichlet(nsamp_per_pop,sample(1:maxscale, replace=T, nclusters)),
                        rdirichlet(nsamp_per_pop,sample(1:maxscale, replace=T, nclusters))), 
                        nrow=(nsamp_per_pop*npop));


simulate_allele_freq <- function(alpha, nSNPs){
  #mu = 1.25e-8
  #alpha = 4*Ne*mu
  beta = alpha
  return(rbeta(nSNPs, shape1 = alpha, shape2 = beta))
}

# we're working with genotype data so equal number of variants across freq range
alpha <- c(1,1,1,1);

# size = nclusters x nSNPs
freq_mat <- t(matrix(unlist(lapply(1:nclusters, function(n) simulate_allele_freq(alpha[n],nSNPs =5))),ncol=nclusters));

simulate_binomial_model <- function(omega, freq_mat)
{
  # nsamp: the number of samples
  # nsites: the number of segregating sites 
  # omega: the topic proportion matrix (nsamp * # of topics)
  # freq_mat: the allele frequencies matrix (# of topics * nsites)
  
  # prod[1,1] = \Sum_{k=1}^{k=nclusters} omega_i,k*f(k, 1)
  # where f(k,1) is the frequency of the 1st SNP for the kth cluster.
  # omega_i,k is the proportion of admixture pop k contributes to individual i
  prod <- omega%*%freq_mat;
  
  sim_data <- apply(prod, c(1,2), function(x) rbinom(1,2,x))
  
  return(sim_data)
}

data <- simulate_binomial_model(omega,freq_mat)
