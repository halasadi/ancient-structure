source('main.R')
library('data.table')
sim_drift<-function(N=40, p0=.2,ngen=100, nSNPs = 1){
  pvec = numeric(nSNPs)
  pvec = p0 	 
  for( i in 2:ngen){
    pvec = rbinom(nSNPs,2*N,pvec)/(2*N)
  }
  return( pvec )
}

sim_ancient_admix <- function(fa1, fa2, t1, t2, alpha, N, nsamp){
  nSNPs <- length(fa1)
  fb1 = sim_drift(p0 = fa1, ngen = t1, N = N, nSNPs = nSNPs)
  fb3 = sim_drift(p0 = fa2, ngen = t1, N = N, nSNPs = nSNPs) 
  fb2 = alpha*fb1 + (1-alpha)*fb3
  
  fc1 = sim_drift(p0 = fb1, ngen = t2, N = N, nSNPs = nSNPs)
  fc2 = sim_drift(p0 = fb2, ngen = t2, N = N, nSNPs = nSNPs)
  fc3 = sim_drift(p0 = fb3, ngen = t2, N = N, nSNPs = nSNPs)
  
  fall = rbind(fc1, fc2, fc3);
  counts_data <- matrix(nrow = 0, ncol = nSNPs)
  for(n in 1:3)
  {
    counts_data = rbind(counts_data,do.call(cbind, lapply(1:nSNPs,function(x) rbinom(nsamp, 2, prob = fall[n,x]))))
  }
  
  return(counts_data)
}

nSNPs = 1000
#fa1 = runif(nSNPs)
#fa2 = runif(nSNPs)
#fa2 = 0.95*fa1 + 0.05*runif(nSNPs, 0, 0.1)

freq_mat <- data.frame(fread("../external_data/ancestral_freqs.txt"), row.names = TRUE);
freq_mat <- freq_mat[1:2,1:nSNPs]
t1 = 10
t2 = 100
alpha = 0.5
nsamp = 100
N = 1e4
data = sim_ancient_admix(freq_mat[1,], freq_mat[2,], t1, t2, alpha, N, nsamp)

K_known = 1
K_unknown = 1
K_pooled = K_known + K_unknown
#freq_mat = rbind(fa1, fa2)
freq_mat = rbind(fa1)
out <- ancient_structure(geno_data = data, f_known = t(freq_mat), K_unknown = K_unknown, max_iter = 500)
barplot(t(out$q),col=2:(K_pooled+1),axisnames=F,space=0,border=NA,main=paste("estd structure: No. of clusters=",K_pooled),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
