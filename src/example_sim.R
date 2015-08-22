
### SQUAREM implementation (Main function)

setwd('/Users/kushal/Documents/ancient-structure/src')
source('simulate_model.R')
source('ancient_structure.R');
## Assignment operations 
npop=5;
nsamp_per_pop=100;
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

# size = nSNPs x nclusters
freq_mat <- matrix(unlist(lapply(1:nclusters, function(n) simulate_allele_freq(alpha[n],nSNPs = 1000))),ncol=nclusters);


data <- simulate_binomial_model(omega,freq_mat)



K_known = 2
K_unknown = 1;
K_pooled <- K_known +K_unknown;
system.time(
out <- ancient_structure(geno_data = data, f_known = freq_mat[,(K_unknown+1):K_pooled], K_unknown = K_unknown, max_iter = 500)
)

barplot(t(omega),col=2:(nclusters+1),axisnames=F,space=0,border=NA,main=paste("true structure: No. of clusters=",nclusters),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
barplot(t(out$q),col=2:(nclusters+1),axisnames=F,space=0,border=NA,main=paste("estd structure: No. of clusters=",nclusters),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
