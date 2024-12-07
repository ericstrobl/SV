posthoc_test <- function(mod,Tx,Y,nperms=10000){
  # post-hoc permutation testing
  # Inputs:
  #   mod = original model from SV function
  #   Tx = vector of treatments for n samples
  #   Y = n by p matrix of p rating scale items
  #   nperms = number of permutations (default: 10000)
  #
  # Outputs:
  #   p = uncorrected p-values for each treatment pair and factor
  #   pFWER = FWER corrected p-values across treatments
  #   stat = difference statistic
  #
  #
  # written by Eric V. Strobl, 12/3/2024
  
  tx = unique(Tx)
  nt = length(tx)
  
  pval = matrix(0,(nt^2-nt)/2,ncol(mod$MR)) # uncorrected p-values
  pFWER = pval
  n = nrow(Y)
  
  # compute difference and absolute difference statistics for original samples
  abs_diff = pval
  diff = pval
  for (j in 1:ncol(mod$MR)){
    abs_diff[,j] = as.numeric(dist(mod$MR[,j]))
    
    z = outer(mod$MR[,j],mod$MR[,j],'-'); 
    diff[,j] = z[lower.tri(z)]
  }
  
  # perform permutation statistic
  for (p in 1:nperms){
    cat('\r',p)
    perm = sample(1:n,n,replace=FALSE)
    mod1 = SV(Tx,Y[perm,],ee = mod$eigen)
    
    range = apply(mod1$MR,2,max) - apply(mod1$MR,2,min) # range statistic for FWER
    for (j in 1:ncol(mod$MR)){
      pFWER[,j] = pFWER[,j] + (abs_diff[,j] < range[j])
      pval[,j] = pval[,j] + (abs_diff[,j] < as.numeric(dist(mod1$MR[,j])))
    }
  }
  pval = pval/nperms
  pFWER = pFWER/nperms
  
  # label comparisons
  labels = c()
  for (c in 1:(length(tx)-1)){
    for (d in (c+1):length(tx)){
      labels = c(labels, paste(as.character(tx[c]), "-", as.character(tx[d]),sep=""))
    }
  }
  rownames(pval) = labels
  rownames(pFWER) = labels
  rownames(diff) = labels
  
  return( list(p = pval, pFWER = pFWER, stat =  diff) )
}
