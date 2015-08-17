
sim_drift <- function(fa1, fa2, t1, t2, alpha, N, nsamp){
  nSNPs <- length(fa1)
  fb1 = rnorm(mean = fa1, sd = sqrt(t1/(2*N)), n = nSNPs)
  fb2 = alpha*fb1 + (1-alpha)*fb3
  fb3 = rnorm(mean = fa2, sd = sqrt(t1/(2*N)), n = nSNPs)
  fc1 = rnorm(mean = fb1, sd = sqrt(t2/(2*N)), n = nSNPs)
  fc2 = rnorm(mean = fb2, sd = sqrt(t2/(2*N)), n = nSNPs)
  fc3 = rnorm(mean = fb3, sd = sqrt(t2/(2*N)), n = nSNPs)
  fall = rbind(fc1, fc2, fc3);
  counts_data <- matrix(nrow = 0, ncol = nSNPs)
  for(n in 1:3)
  {
    counts_data = rbind(counts_data,do.call(cbind, lapply(1:nSNPs,function(x) rbinom(nsamp, 2, prob = fall[n,x]))))
  }
  
}