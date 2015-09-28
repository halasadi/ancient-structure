
### SQUAREM implementation (Main function)
set.seed(10)
source('simulate_model.R')
source('ancient_structure.R');
## Assignment operations
npop=3;
nsamp_per_pop=100;
pop_labs <- unlist(lapply(1:npop, function(n) rep(toString(n),nsamp_per_pop)));
source_labs <- toString(1);
nclusters = 3
maxscale = 10
# omega is defined to be the admixture proportions
# size = nsamp * nclusters

omega <- matrix(data=0,nrow=0,ncol=nclusters);
for(pop in 1:npop) {
omega <- rbind(omega,rdirichlet(nsamp_per_pop,sample(1:maxscale, replace=T, nclusters)));
}

omega <- fix_clus_mem(omega, pop_labs, source_labs);

simulate_allele_freq <- function(alpha, nSNPs){
  #mu = 1.25e-8
  #alpha = 4*Ne*mu
  beta = alpha
  return(rbeta(nSNPs, shape1 = alpha, shape2 = beta))
}

# we're working with genotype data so equal number of variants across freq range
alpha <- c(1,1,1,1);
nSNPs <- 1000;
# size = nSNPs x nclusters
freq_mat <- matrix(unlist(lapply(1:nclusters, function(n) simulate_allele_freq(alpha[n],nSNPs = nSNPs))),ncol=nclusters);


data <- simulate_binomial_model(omega,freq_mat)



K_known = 1;
K_unknown = 2;
K_pooled <- K_known +K_unknown;

f_obs <- matrix(0,nrow=nSNPs, ncol=K_known);
n_obs = rep(2,K_known);
if(min(n_obs)!=0){
for (k in 1:K_known){
  f_obs[,k] = rbinom(nSNPs,2*n_obs[k], freq_mat[,k])/(2*n_obs[k])
}
}

system.time(
out <- ancient_structure(geno_data = data, f_obs = f_obs, n_obs = n_obs, 
                         pop_labs=pop_labs, source_labs=source_labs,
                         K_unknown = K_unknown, max_iter = 600)
)

barplot(t(omega),col=2:(nclusters+1),axisnames=F,space=0,border=NA,main=paste("true structure: No. of clusters=",nclusters),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

docweights <- out$q;

library(permute);
library("BioPhysConnectoR");
perm_set=rbind(1:nclusters,allPerms(1:nclusters));
diff=array(0,dim(perm_set)[1]);
for (p in 1:dim(perm_set)[1])
{
  temp=docweights[,perm_set[p,]];
  diff[p]=fnorm(temp,omega);
}

p_star=which(diff==min(diff));
docweights=docweights[,perm_set[p_star,]];


barplot(t(docweights),col=2:(nclusters+1),axisnames=F,space=0,border=NA,main=paste("estd structure: No. of clusters=",nclusters),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)
