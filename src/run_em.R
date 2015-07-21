
## Topic model fit with known ancestral allele frequencies and unknown ancestral allele 
## frequencies 

update <- function(q_in, f_unknown_in, f_known_in, geno_data_local)
  # q_in: topic proportions (n x k)
  # f_unknown: allele frequences of the unknown ancestral populations (nSNPs x k_unknown)
  # f_knowm: allele frequencies of the known ancestral populations (nSNPs x k_known)
  # geno_data_local <- the genotype matrix data (nsamp X nSNPs)
  # if k_unknown =0, use the input for f_unknown to be matrix(nrow=nSNPs,ncol=0)
  
{
  f_pooled_in <- cbind(f_unknown_in, f_known_in);
  nSNPs = dim(f_pooled_in)[1]
  nsamp = dim(q_in)[1]
  K_pooled = dim(f_pooled_in)[2] 
  
  K_known <- dim(f_known_in)[2];
  K_unknown <- dim(f_unknown_in)[2];
  
  a <- array(0, c(nsamp, nSNPs, K_pooled));
  
  for(i in 1:nsamp)
  {
    for(j in 1:nSNPs)
    {
      a[i,j,] = geno_data_local[i,j]* (q_in[i,]*f_pooled_in[j,])/(q_in[i,]%*%f_pooled_in[j,]);
    }
  }
  
  
#  a_outer<- array(0, c(dim(q_in)[1], dim(f_pooled_in)[1], dim(f_pooled_in)[2]));
  
#  for(i in 1:dim(q_in)[1])
#  {
#    a_outer[i,,] = (geno_data_local[i,]*(q_in[i,]*f_pooled_in))/matrix(rep(t(q_in[i,]%*%t(f_pooled_in)),dim(f_pooled_in)[2]),nrow=dim(f_pooled_in)[1])
#  }
  
  b <- array(0, c(nsamp, nSNPs, K_pooled));
  
  for(i in 1:nsamp)
  {
    for(j in 1:nSNPs)
    {
      b[i,j,] = (2-geno_data_local[i,j])* (q_in[i,]*(1-f_pooled_in[j,]))/(q_in[i,]%*%(1-f_pooled_in[j,]));
    }
  }
  
  if (K_unknown != 0){
    f_unknown_out <- array(0, c(nSNPs,K_unknown));
    
    for (j in 1:nSNPs)
    {
      for(k in 1:K_unknown){
        f_unknown_out[j,k] <- sum(a[,j,k])/(sum(a[,j,k])+sum(b[,j,k]));
      }
    }
    
    f_known_out <- f_known_in; 
  }
  
  q_out <- array(0, c(nsamp,K_pooled))
  
  for (i in 1:nsamp)
  {
    for(k in 1:K_pooled){
      q_out[i,k] <- 0.5 * mean(a[i,,k]) + 0.5* mean(b[i,,k]);
    }
  }
  
  
  out <- list("f_unknown"=f_unknown_out,"f_known"=f_known_out,"q"=q_out);
  return(out)
  
}

topic_model_ancestral <- function(geno_data, freq_known, n_unknown)
{
  
  
  
}