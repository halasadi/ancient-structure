
# @: function to simulate a N \times L matrix for L sites and N samples 
# day populations and N_a ancient populations. 

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

