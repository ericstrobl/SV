# Introduction

This is the R code repository for the Supervised Varimax (SV) algorithm. Treatments tend to diffusely affect many items of a clinical rating scale with large amounts of noise. SV thus takes the individual items of a rating scale and compresses them into a few (e.g., 2-5) factors that maximally differentiate between the treatments using randomized clinical trial data. The algorithm specifically performs a Varimax rotation on the matrix of causal effects from treatment to factors so that the causal effects are sparse and different.

# Installation

> library(devtools)

> install_github("ericstrobl/SV")

> library(SV)

# Inputs

`Tx` = vector of treatments for n samples

`Y` = n by p matrix of p rating scale items

`nc` (optional) = number of components/factors; default is number of unique treatments

`eigen` (optional) = eigendecomposition of cor(Y) from previous output of SV; useful for minimizing run-time of permutation tests; default is NULL

# Run the Algorithm

> data = generate_synth(nsamps=1000, nF=4) # generate synthetic data with five treatments and four latent factors

> out=SV(Tx=data$Tx,Y=data$Y) # run the SV algorithm on the data, where mod$Tx is a vector of treatments, and mod$Y are the individual items of a rating scale

# Outputs

A list with:

`MR` = causal effects from treatments to factors

`RtW` = causal effects from factors to treatments

`tx_effects` = (MR) %*% (RtW), or the causal effects from treatments to factors

`optimal_outcomes` = rotated factors

`R` = the rotation matrix found by Varimax

`sgn` = sign flip done to ensure sign determinancy

`io` = re-ordering of optimal outcomes done to ensure permutation determinancy

`eigen` = results of eigendecomposition of cor(Y)

# Diagnostics

> plot(sort(apply(SV$loadings^2,2,var),decreasing=T)) # cut at the elbow point, should be 4 factors in this example

# Permutation tests

Omnibus:

> aa

Post-hoc:

> aa
