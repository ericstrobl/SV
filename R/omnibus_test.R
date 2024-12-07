omnibus_test <- function(mod,Tx,Y,nperms=10000){
  # omnibus permutation testing
  # Inputs:
  #   mod = original model from SV function
  #   Tx = vector of treatments for n samples
  #   Y = n by p matrix of p rating scale items
  #   nperms = number of permutations (default: 10000)
  #
  # Outputs:
  #   p = p-value
  #   stat = absolute sum statistic
  #
  #
  # written by Eric V. Strobl, 12/3/2024

  pval = 0
  stat = sum(abs(mod$MR))
  n = nrow(Y)
  for (p in 1:nperms){
    cat('\r',p)
    perm = sample(1:n,n,replace=FALSE)
    mod1 = SV(Tx,Y[perm,],ee = mod$eigen)
    pval = pval + (stat <= sum(abs(mod1$MR)))
  }

  return( list(p = pval/nperms, stat =  stat) )
}
