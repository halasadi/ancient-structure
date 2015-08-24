
###  Transformation and Reverse transform functions


## the reverse transform function takes a simplex vector and un-simplexes it on (-infty,infty)

reverse_transform=function(x) 
{
  return(log((abs(x[2:length(x)])+1e-7)/(abs(x[1])+1e-7)));
}

# the transform function simplexes a vector

transform <- function(y) 
{
  temp =c(1,exp(y));
  out=temp/sum(temp);
  return(out)
}

### loglik_squarem : performs log-likelihood of the data under Binomial model of admixture

# @param_vec_in : a vector of size nsamp*n_clus+n_clus_known*n_genes+n_clus_unknown*ngenes: rev_trans(q), logit(f_known), logit(f_unknown) vectorized

loglik_squarem <- function(param_vec_in, geno_data, nsamp, K_unknown, K_known, nSNPs)
{
  K_pooled <- K_known + K_unknown;
  
  # rev_q_in is nsamp*(K_pooled-1) matrix, each row is a reverse transform on q (topic prop matrix)
  # so has one less co-ordinate as it is not constrained
  
  rev_q_in = matrix(param_vec_in[(1:(nsamp*(K_pooled-1)))],nrow = nsamp, ncol = (K_pooled-1));
  
  
  # transforming the rev_q_in to q_in by adding the sum=1 constraint
  
  q_in <- t(apply(rev_q_in, 1,function(x) transform(x)));
  
  
  # extracting the unknown allele frequencies data from the param_vec_in object
  
  temp <- param_vec_in[-(1:(nsamp*(K_pooled-1)))];
  f_unknown_in <- inv.logit(matrix(temp[0:(nSNPs*K_unknown)], nrow=nSNPs, ncol=K_unknown))
  
  
  # extracting the known allele frequencies data from the param_vec_in object
  
  beg = nSNPs*K_unknown + 1;
  end = nSNPs*K_pooled;
  f_known_in <- inv.logit(matrix(temp[beg:end], nrow=nSNPs, ncol=K_known));
  
  # pooling the unknown and known allele frequencies
  
  f_pooled_in <- cbind(f_unknown_in, f_known_in);
  
  # loglikelihood computation
  prod <- q_in %*% t(f_pooled_in);
  return( -sum(log(prod)*geno_data - (2-geno_data)*log(1-prod)))
}

library(SQUAREM)


# SQUAREM wrapper function

# @param_vec_in : a vector of size nsamp*n_clus+n_clus_known*n_genes+n_clus_unknown*ngenes: rev_trans(q), logit(f_known), logit(f_unknown) vectorized
# update_EM is the main function 
# This function prepares the input for update_EM and then reverse transforms the output for squarem interpolation



update_squarem <- function(param_vec_in, geno_data, nsamp, K_unknown, K_known, nSNPs)
{
  
  K_pooled <- K_known +K_unknown;
  
  # rev_q_in is nsamp*(K_pooled-1) matrix, each row is a reverse transform on q (topic prop matrix)
  # so has one less co-ordinate as it is not constrained
  
  rev_q_in = matrix(param_vec_in[(1:(nsamp*(K_pooled-1)))],nrow = nsamp, ncol = (K_pooled-1));
  
  # transforming the rev_q_in to q_in by adding the sum=1 constraint
  
  q_in <- t(apply(rev_q_in, 1,function(x) transform(x)));
  
  # extracting the unknown allele frequencies data from the param_vec_in object and inv transforming it
  
  temp <- param_vec_in[-(1:(nsamp*(K_pooled-1)))];
  f_unknown_in <- inv.logit(matrix(temp[0:(nSNPs*K_unknown)], nrow=nSNPs, ncol=K_unknown));
  
  # extracting the known allele frequencies data from the param_vec_in object and inv transforming it
  
  beg = nSNPs*K_unknown + 1;
  end = nSNPs*K_pooled;
  f_known_in <- inv.logit(matrix(temp[beg:end], nrow=nSNPs, ncol=K_known));
  
  # using the EM update scheme- the main function 
  out <- update_EM(q_in, f_unknown_in,f_known_in,geno_data);
  rev_q_out <- as.matrix(t(apply(out$q, 1, function(x) reverse_transform(x))));
  param_vec_out <- c(as.vector(rev_q_out),as.vector(logit(out$f_unknown)),as.vector(logit(out$f_known)));
  return(param_vec_out)
}

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
  
  
# 
#  a_outer <- mclapply(1:nsamp, function(i) sweep((geno_data[i,]*(matrix(rep(q_in[i,],nSNPs),nrow=nSNPs,byrow=TRUE)*f_pooled_in)),1,(q_in[i,]%*%t(f_pooled_in)),'/'));
#  for(n in 1:nsamp)
#  {
#    a[n,,] <- a_outer[[n]];
#  }
  
  a <- array(0, c(nsamp, nSNPs, K_pooled));
  
  a_outer <- mclapply(1:nsamp, function(i) apply((geno_data[i,]*(matrix(rep(q_in[i,],nSNPs),nrow=nSNPs,byrow=TRUE)*f_pooled_in)),2,'/',(q_in[i,]%*%t(f_pooled_in))),mc.cores = detectCores());
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
  
  b_outer <- mclapply(1:nsamp, function(i) apply(((2-geno_data[i,])*(matrix(rep(q_in[i,],nSNPs),nrow=nSNPs,byrow=TRUE)*(1-f_pooled_in))),2,'/',(q_in[i,]%*%t((1-f_pooled_in)))),mc.cores=detectCores());
  
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


