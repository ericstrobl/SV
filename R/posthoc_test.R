posthoc_test <- function(mod,Tx,Y,nE,nperms=10000){
  # post-hoc permutation testing
  #   mod = original model from SV function
  #   Tx = vector of treatments for n samples
  #   Y = n by p matrix of p rating scale items
  #   nE = number of components/factors ELIMINATED after diagnostics
  #
  #
  # written by Eric V. Strobl, 12/3/2024
  
  require(qvalue)
  
  pval = 0
  stat = (apply(mod$MR,2,max) - apply(mod$MR,2,min))
  n = nrow(Y)
  for (p in 1:nperms){
    cat('\r',p)
    perm = sample(1:n,n,replace=FALSE)
    mod1 = SV(Tx,Y[perm,],ee = mod$eigen)
    pval = pval + (stat <= (apply(mod1$MR,2,max) - apply(mod1$MR,2,min)))
  }
  
  ps = pval/nperms
  d = length(ps)
  qs = qvalue(ps,pi0=(d-nE-1)/d)$qvalues
  
  return( list(p = pval/nperms, p_FDR = qs, stat =  stat) )
}
