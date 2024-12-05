
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

nc = ncol(Fx)
coef = sweep(T_F %*% F_Y,2,mY,FUN="-") %*% diag(1/sdY)
# print(coef)

ee = eigen(cov(Y))

# print(ee$values[1:nc])
R = my_varimax(coef%*% ee$vectors[,1:nc,drop=FALSE] %*% diag(1/sqrt(ee$values[1:nc])),normalize=FALSE)$rotmat
# R = diag(nc) # compare against normal supervised PCA
factors = Y %*% ee$vectors[,1:nc,drop=FALSE] %*% diag(1/sqrt(ee$values[1:nc])) %*% R

sgn = c(sign(cor(factors,rowSums(Y)))) # eliminate sign indeterminancy
io = order(cor(factors,rowSums(Y))^2,decreasing=T) # eliminate permutation indeterminancy, order by amount of variance explained

factors = sweep(factors,2,sgn,FUN="*")[,io]

# compute loadings (i.e., tx effects on factors)
loadings = coef%*% ee$vectors[,1:nc,drop=FALSE] %*% diag(1/(sqrt(ee$values[1:nc]))) %*% R
rownames(loadings) = txs
loadings = sweep(loadings,2,sgn,FUN="*")[,io]

# compute factor effects on Y
factor_effects = t(R) %*% diag(sqrt(ee$values[1:nc])) %*% t(ee$vectors[,1:nc])
factor_effects = sweep(factor_effects,1,sgn,FUN="*")[io,]

# compute tx effects on Y
tx_effects = loadings %*% factor_effects

return(list(Y=Y, loadings = loadings, eigen = ee, factors = factors, rotat = R,
            factor_effects = factor_effects,
            tx_effects = tx_effects, Tx = Tx,vari_input = coef%*% ee$vectors[,1:nc,drop=FALSE] %*% diag(1/sqrt(ee$values[1:nc]))))

}
