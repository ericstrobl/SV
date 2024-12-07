# Introduction

This is the R code repository for the Supervised Varimax (SV) algorithm. Treatments tend to diffusely affect many items of a clinical rating scale with large amounts of noise. SV thus takes the individual items of a rating scale and compresses them into a few factors that maximally differentiate between the treatments using randomized clinical trial data. The algorithm specifically performs a Varimax rotation on the matrix of causal effects from treatment to factors so that the causal effects are sparse and different.

The academic paper describing SV in detail can be found [here](https://www.medrxiv.org/content/10.1101/2024.12.03.24318424v1); it is currently under review. All code was tested in R version 4.3.1.

# Installation
First install BiocManager. Then:

> library(devtools)

> install_github("ericstrobl/SV")

> library(SV)

# Inputs

`Tx` = vector of treatments for n samples

`Y` = n by p matrix of p rating scale items

`nc` (optional) = number of components/factors; default is number of unique treatments

`ee` (optional) = eigendecomposition of cor(Y) from previous output of SV; useful for minimizing run-time of permutation tests; default is NULL

# Run the Algorithm

> data = generate_synth(nsamps=10000, nF=3) # generate synthetic data with five treatments and three latent factors

> mod=SV(Tx=data$Tx,Y=data$Y) # run the SV algorithm on the data, where data$Tx is a vector of treatments, and data$Y are the individual items of a rating scale

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

> plot(sort(apply(mod$MR^2,2,var),decreasing=T)) # cut at and below the elbow point; the elbow point is very close to zero and should be at the 4th factor in this case, so there are 3 clinically significant factors

# Permutation tests

Omnibus:

> omnibus_test(mod,data$Tx,data$Y)

Post-hoc:

> posthoc_test(mod,data$Tx,data$Y)
