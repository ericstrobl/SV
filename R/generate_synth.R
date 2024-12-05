
generate_synth <- function(nF = sample(2:4,1),nTx = 5, nY = 20, nsamps=1000){

Tx = rep(1:5,length.out=nsamps)
txs = unique(Tx)
t_idx = matrix(0,nsamps,length(txs))
for (t in 1:length(txs)){
  t_idx[which(Tx == txs[t]),t] = 1
}

T_F = matrix((runif(nTx * nF)*0.75+0.25)*sample(c(1,-1),nTx*nF,replace=TRUE),nTx,nF)
F_Y = matrix((runif(nF * nY)*0.75+0.25)*sample(c(1,-1),nF*nY,replace=TRUE),nF,nY)
E_f = matrix(rt(nsamps*nF,df=3),nsamps,nF)

Fx = t_idx %*% T_F + E_f
Y = Fx %*% F_Y
mY = apply(Y,2,mean)
sdY = apply(Y,2,sd)
Y = sweep(Y,2,1/sdY,FUN='*')

return(list(Y=Y, Tx = Tx))

}
