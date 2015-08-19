
## Topic model fit with known ancestral allele frequencies and unknown ancestral allele 
## frequencies 

update_EM <- function(q_in, f_unknown_in, f_known_in, geno_data)
  # q_in: topic proportions (n x k)
  # f_unknown: allele frequences of the unknown ancestral populations (nSNPs x k_unknown)
  # f_knowm: allele frequencies of the known ancestral populations (nSNPs x k_known)
  # geno_data <- the genotype matrix data (nsamp X nSNPs)
  # if k_unknown =0, use the input for f_unknown to be matrix(nrow=nSNPs,ncol=0)
  
{
  f_pooled_in <- cbind(f_unknown_in, f_known_in);
  nSNPs = dim(f_pooled_in)[1]
  nsamp = dim(q_in)[1]
  K_pooled = dim(f_pooled_in)[2] 
  
  K_known <- dim(f_known_in)[2];
  K_unknown <- dim(f_unknown_in)[2];
  
###   Deriving a using for loop without vectorization and paralleization (commented)
#  a_ser <- array(0, c(nsamp, nSNPs, K_pooled));
  
#  for(i in 1:nsamp)
#  {
#    for(j in 1:nSNPs)
#    {
#      a_ser[i,j,] = geno_data[i,j]* (q_in[i,]*f_pooled_in[j,])/(q_in[i,]%*%f_pooled_in[j,]);
#    }
#  }
  
###   Deriving a using  vectorization and paralleization 
  
  
  a <- array(0, c(nsamp, nSNPs, K_pooled));
  a_outer <- mclapply(1:nsamp, function(i) sweep((geno_data[i,]*(matrix(rep(q_in[i,],nSNPs),nrow=nSNPs,byrow=TRUE)*f_pooled_in)),1,(q_in[i,]%*%t(f_pooled_in)),'/'));
  for(n in 1:nsamp)
  {
    a[n,,] <- a_outer[[n]];
  }
  
  a_outer <- mclapply(1:nsamp, function(i) apply((geno_data[i,]*(matrix(rep(q_in[i,],nSNPs),nrow=nSNPs,byrow=TRUE)*f_pooled_in)),2,'/',(q_in[i,]%*%t(f_pooled_in))),mc.cores=detectCores());
  for(n in 1:nsamp)
  {
    a[n,,] <- a_outer[[n]];
  }
###   Deriving b using for loop without vectorization and paralleization (commented)
  
#  b_ser <- array(0, c(nsamp, nSNPs, K_pooled));
  
#  for(i in 1:nsamp)
#  {
#    for(j in 1:nSNPs)
#    {
#      b_ser[i,j,] = (2-geno_data[i,j])* (q_in[i,]*(1-f_pooled_in[j,]))/(q_in[i,]%*%(1-f_pooled_in[j,]));
#    }
#  }
  
  
## Deriving b using vectorization on SNPs and parallelization on the samples 
  
  b <- array(0, c(nsamp, nSNPs, K_pooled));
  
  b_outer <- mclapply(1:nsamp, function(i) apply(((2-geno_data[i,])*(matrix(rep(q_in[i,],nSNPs),nrow=nSNPs,byrow=TRUE)*(1-f_pooled_in))),2,'/',(q_in[i,]%*%t((1-f_pooled_in)))), mc.cores=detectCores());
  
  for(n in 1:nsamp)
  {
    b[n,,] <- b_outer[[n]];
  }
  
  f_unknown_out <- matrix(nrow=nSNPs,ncol=0)
  f_known_out <- f_known_in;
  
  if (K_unknown != 0){
    f_unknown_out <- sapply(1:K_unknown, function(k) colSums(a[,,k])/(colSums(a[,,k])+colSums(b[,,k])))
#    f_unknown_out <- array(0, c(nSNPs,K_unknown));
#    for (j in 1:nSNPs)
#    {
#      for(k in 1:K_unknown){
#        f_unknown_out[j,k] <- sum(a[,j,k])/(sum(a[,j,k])+sum(b[,j,k]));
#      }
#    }
    
  }
  
  
  
  q_out <- t(sapply(1:nsamp, function(i) 0.5 * colMeans(a[i,,]) + 0.5* colMeans(b[i,,])));
  out <- list("f_unknown"=f_unknown_out,"f_known"=f_known_out,"q"=q_out);
  return(out)
  
}

