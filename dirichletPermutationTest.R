#######################################
## This code implements the Diriclet
##  likelihood ratio change point test
##  in Prabuchandran (2019)
##
## There, we search over all possible change points
##  in a sequence of iid Dirchlet data. The shape
##  parameter is estimated under the null (no change) and
##  alternative (change occurred, two shape terms).
##  The difference in log likelihoods is calculated and the maximum
##  difference is the test statistic.
## The sampling distribution of the statistic is found through permutation
##  methods, we perform 1000 permutations here.
##
## We also implement a modified version of the test where we allow
##  the data to be seasonal. We find a shape parameters within seasons
##  and permutate within season.
## 
## You will also find a multiple change point version of the test implemented
##  using the algorithm in Prabuchandran (2019).
##
## Find an estimate for the shape can be difficult in small samples
##  thus we needed to make some limitation on locations of potential
##  change points. We also needed to do some exception handling in 
##  R with the tryCatch() function.

library(DirichletReg)
library(sirt)
library(tidyverse)


parametricDirichletTestSeasons <- function(data) {
  data <- as.matrix(data)
  seasons <- rep(1:3, ceiling(dim(data)[1]/3))[1:dim(data)[1] ]
  spring <- which(seasons==1)
  summer <- which(seasons==2)
  fall <- which(seasons==3)

  null_llik <- function(data) {
    spr.alpha <- dirichlet.mle(data[spring,])$alpha
    sum.alpha <- dirichlet.mle(data[summer,])$alpha
    fal.alpha <- dirichlet.mle(data[fall,])$alpha
    alpha.mat <- matrix(nrow=dim(data)[1], ncol=dim(data)[2], byrow=TRUE, 0)
    alpha.mat[spring,] <- matrix(spr.alpha, ncol=4, nrow=length(spring), byrow=TRUE)
    alpha.mat[summer,] <- matrix(sum.alpha, ncol=4, nrow=length(summer), byrow=TRUE)
    alpha.mat[fall,] <- matrix(fal.alpha, ncol=4, nrow=length(fall), byrow=TRUE)
    ddirichlet(as.matrix(data), alpha=alpha.mat, log=TRUE, sum.up=TRUE)
  }
  
  alt_llik <- function(data, tau) {
    spr.alpha1 <- dirichlet.mle(data[spring[spring<=tau],])$alpha
    sum.alpha1 <- dirichlet.mle(data[summer[summer<=tau],])$alpha
    fal.alpha1 <- dirichlet.mle(data[fall[fall<=tau],])$alpha
    spr.alpha2 <- dirichlet.mle(data[spring[spring>tau],])$alpha
    sum.alpha2 <- dirichlet.mle(data[summer[summer>tau],])$alpha
    fal.alpha2 <- dirichlet.mle(data[fall[fall>tau],])$alpha
    
    alpha.mat <- matrix(nrow=dim(data)[1], ncol=dim(data)[2])
    alpha.mat[spring[spring<=tau],] <- matrix(spr.alpha1, ncol=4, nrow=sum(spring<=tau), byrow=TRUE)
    alpha.mat[spring[spring >tau],] <- matrix(spr.alpha2, ncol=4, nrow=sum(spring>tau), byrow=TRUE)

    alpha.mat[summer[summer<=tau],] <- matrix(sum.alpha1, ncol=4, nrow=sum(summer<=tau), byrow=TRUE)
    alpha.mat[summer[summer >tau],] <- matrix(sum.alpha2, ncol=4, nrow=sum(summer>tau), byrow=TRUE)

    alpha.mat[fall[fall<=tau],] <- matrix(fal.alpha1, ncol=4, nrow=sum(fall<=tau), byrow=TRUE)
    alpha.mat[fall[fall >tau],] <- matrix(fal.alpha2, ncol=4, nrow=sum(fall>tau), byrow=TRUE)

    ddirichlet(as.matrix(data), alpha=alpha.mat, log=TRUE, sum.up=TRUE)
  }
  
  ## We force at least 3 years in a given window so we can get
  ## MLE estimates. 
  out <- sapply(13:(dim(data)[1]-12), alt_llik, data=data)
  tau_hat <- (13:(dim(data)[1]-12))[which.max(out)]   ## Chosen ch-pt from data
  z_star <- max(out) - null_llik(data)              ## Z-star, test stat from data
  
  get_sample <- function(tau) {
    perm_data <- matrix(nrow=dim(data)[1], ncol=dim(data)[2], 0)
    
    ## Shuffle within each season
    ind <- sample(spring, size=length(spring), replace=FALSE)
    perm_data[spring,] <- data[ind,]
    ind <- sample(summer, size=length(summer), replace=FALSE)
    perm_data[summer,] <- data[ind,]
    ind <- sample(fall, size=length(fall), replace=FALSE)
    perm_data[fall,] <- data[ind,]
    
    ### Here we search over all possible change points within shuffled data
    ###   Find the likelihood value and pick the largest
    ### We return the difference in log-likelihoods, effectively
    ###   a sample of z-star
    out <- sapply(13:(dim(data)[1]-12), alt_llik, data=perm_data)
    max(out) - null_llik(perm_data)
  }

  ### Occasionally get a singularity error in the log likelihood calculation
  ### due to a poor shuffled sample, so a little error handling.
  ###
  ### This method effectively "cheats" a bit, because it uses the picked
  ###   tau-hat value, and you use it based on the permutations
  z_star_sample <- sapply(1:1100, function(x) {
    tryCatch({get_sample(tau_hat)}, error=function(err){return(NA)}) } )
  z_star_sample <- sample(na.omit(z_star_sample), size=1000, replace=FALSE)
  c(tau_hat, z_star, (sum(z_star_sample>z_star)+1)/1001)
}



parametricDirichletTest <- function(data) {
  data <- as.matrix(data)
  ind <- 1:dim(data)[1]

  null_llik <- function(data) {
    alpha.hat <- dirichlet.mle(data)$alpha
    ddirichlet(as.matrix(data), alpha=alpha.hat, log=TRUE, sum.up=TRUE)
  }
  
  alt_llik <- function(data, tau) {
    alpha1 <- dirichlet.mle(data[ind[ind<=tau],])$alpha
    alpha2 <- dirichlet.mle(data[ind[ind >tau],])$alpha
   
    alpha.mat <- matrix(nrow=dim(data)[1], ncol=dim(data)[2])
    alpha.mat[ind[ind<=tau],] <- matrix(alpha1, ncol=4, nrow=sum(ind<=tau), byrow=TRUE)
    alpha.mat[ind[ind >tau],] <- matrix(alpha2, ncol=4, nrow=sum(ind >tau), byrow=TRUE)
    
    ddirichlet(as.matrix(data), alpha=alpha.mat, log=TRUE, sum.up=TRUE)
  }
  
  ## We force at least 3 observations in a given window so we can get
  ## MLE estimates. 
  out <- sapply(5:(dim(data)[1]-4), alt_llik, data=data)
  tau_hat <- (5:(dim(data)[1]-4))[which.max(out)]   ## Chosen ch-pt from data
  z_star <- max(out) - null_llik(data)              ## Z-star, test stat from data
  
  get_sample <- function(tau) {
    perm_data <- matrix(nrow=dim(data)[1], ncol=dim(data)[2], 0)
    
    ## Shuffle within each season
    ind2 <- sample(ind, size=length(ind), replace=FALSE)
    perm_data[ind,] <- data[ind2,]

    ### Here we search over all possible change points within shuffled data
    ###   Find the likelihood value and pick the largest
    ### We return the difference in log-likelihoods, effectively
    ###   a sample of z-star
    out <- sapply(5:(dim(data)[1]-4), alt_llik, data=perm_data)
    max(out) - null_llik(perm_data)
  }
  
  ### Occasionally get a singularity error in the log likelihood calculation
  ### due to a poor shuffled sample, so a little error handling.
  ###
  ### This method effectively "cheats" a bit, because it uses the picked
  ###   tau-hat value, and you use it based on the permutations
  z_star_sample <- sapply(1:1100, function(x) {
    tryCatch({get_sample(tau_hat)}, error=function(err){return(NA)}) } )
  z_star_sample <- sample(na.omit(z_star_sample), size=1000, replace=FALSE)
  c(tau_hat, z_star, (sum(z_star_sample>z_star)+1)/1001)
}




parametricDirichletMultCP <- function(data) {
  ##
  ## Note, in the permutation test, we need at least
  ## 3 observations for each season in each regime
  ## so we need to be careful with how we implement this
  
  n <- dim(data)[1]
  L0 <- 30    ## This privides an initial search window of 9 to 21.. for a change point
  L <- L0
  batch <- 9  ## Want L0 and batch divisible by 3
  begin <- 1
  cpt <- NULL
  run <- TRUE
  while(run) {
    if(L >= n) {
      L <- n
      run=FALSE
    }
    ## Need at least 21 observations in the window to even try
    ##   We implement this because in small samples calculating
    ##   the likelihood value will fail on occassion.
    if(L-begin>=24) 
      tmp <- parametricDirichletTest(data[begin:L,])
    else tmp<-c(0,0,1)
    if(tmp[3]<0.05) {
      cpt <- c(cpt, tmp[1] )
      begin <- L+1
      L <- L+L0
      
    } else {
      L <- L + batch
    }
  }
  cpt
}


