##
## This code fits the two-state HMM
##   where we separate model the location & scale
##   of the Dirichlet distribution
##
## This code has been streamlined for readability
##   and to demonstrate how to fit the HMM
## 
## The example here is with the phytoplankton data
##   and should correspond to the results presented
##   in the primary manuscript.

#######################
## Preliminary stuff
iter=100000
thin=100
seed=314159

set.seed(seed)

library(tidyverse)
library(rstan)
library(parallel)

rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores(), width=150)

stan_mod_loc_scale <- stan_model(file="../stan_models/dir_reg_loc_and_scale_hmm.stan")

plankton.data <- read.csv("../data_processing/biomass3.csv")
Y <- plankton.data %>%
  dplyr::select(Prop.BG, Prop.Green, Prop.Flag, Prop.Diatom)

#################################
## Setup the stan data
##   A list with all necessary
##   parameters
plankton_stan_data <- within(list(), {
  T <- dim(Y)[1]
  M <- 2
  K <- dim(Y)[2]
  P <- 3
  y <- as.matrix(Y)
  X <- as.matrix(model.matrix(~as.factor(rep(1:3, T/3))) )
  rho_a <- 9.5
  rho_b <- 0.5
  b_loc_sig <- 2.0
  b_scale_sig <- 2.0
})

###################
## initial values for the MCMC algorithm
## 
## list of length 2 since two chains
init.vals.loc.scale <- list()
for(i in 1:2) {
  init.vals.loc.scale[[i]] <- list(prob_remain = array(0.9, dim=c(1)),
                                   b_loc=array(rep(0, 24), dim=c(2,3,4)),
                                   b_scale=array(rep(0,6), dim=c(2,3) ) )
}

#################
## Fetch a posterior sample
fitHMM_loc_and_scale_plankton <- sampling(stan_mod_loc_scale, plankton_stan_data, iter=iter, thin=thin, chains=2, init=init.vals.loc.scale,
                                          cores=2, verbose=TRUE, seed=seed, save_warmup=FALSE) 

###############
## Save rstan result -- pretty big file
##  We extract() values for plots & results elsewhere
save( 
     fitHMM_loc_and_scale_plankton,
     file="plankton_rstan_fits_2state.RData")

