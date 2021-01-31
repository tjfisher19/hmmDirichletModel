##
## This code generates 200x4 simulated datasets
##   from various Dirichlet distributions
##
## We follow the convention that the shape parameter (alpha) of the
##   Dirichlet distribtuion can be decomposed in a location part (theta)
##   and a spread part (tau), thus alpha = tau*theta where tau is a scalar
##   and theta ia p-dimensional vector on the simplex.
## 
## The most basic is essentially a 4-dimensional Uniform-type distribution
## 
## We also generate a dataset with no change in parameters following
##   our plankton data values. That is, we estimate a shape vector
##   for spring, summer and fall, and generate a seasonal Dirichlett
##   dataset with no change in parameters
##
## A seasonal Dirichlett process is generated with a shift in parameters
##   at time point 31 is generated. Up to time point 30, the location parameters follow
##   the first 7 years of plankton dataset (fit with MLE). After point 30, the
##  location parameters are generated using the estimates from the last 7 years.
##
## A dataset is generated with a shift in the spread parameter
##   Here, we follow the same type of setup as above for the location shift
##   except the location parameters stay constant in time (except for seasonal variation)
##   but the spread parameter changes at time point 31.
## 
##   We explored two versions of this change
##     One was based on the first and last 7 year estimates from our real data
##       We found this did not result in much change and all methods struggled
##     The second approach +/- the MLE values for the spread of the full dataset
##       This allowed us to see a difference in the methods and demonstrated
##       the functionality that was proposed
##
## We also generate data where two change points are present (one in location and
##   one in the scale term).
##
## In a covariate example we generate a seasonal time series and use it as
##   the covariate to generate a Dirichlet response. The coefficients
##   are based on the fitted covariate model in the paper and a change point
##   occurs in the relationship
##
## And lastly we generate data where both the location and scale change, but 
##   in such a way that the overall shape does not change too much, thus
##   it may be hard to detect the shift.
##
## Some notes about code, otherwise fairly straightward.
##
##  gtools is used for rdirichlet
##  sirt is used to find MLEs of Dirichlet data
##  tidyverse for data processing

library(gtools)
library(sirt)
library(tidyverse)


## Generate some 4-dimensional uniform type data
## No frills here
uniformNoChangeDataFun <- function(n=63) {
  Y.sim <- rdirichlet(n, alpha=rep(1,4))
  Y.sim
}
set.seed(314159)
uniformNoChangeData <- lapply(rep(63, 200), uniformNoChangeDataFun)


## The remaining of our generated datasets will be based
## on the values estimated from our real data
## We use our plankton data as the baseline
plankton.data <- read.csv("biomass3.csv")
Y <- plankton.data %>%
  dplyr::select(Prop.BG, Prop.Green, Prop.Flag, Prop.Diatom)

####
## First we get the MLE estimates for Spring, Summer and Fall
## Under the assumption of no shifts, just seasonal components
sprg <- seq(1,61,3)
summ <- seq(2,62,3)
fall <- seq(3,63,3)
fit.sprg <- dirichlet.mle(Y[sprg,])
fit.summ <- dirichlet.mle(Y[summ,])
fit.fall <- dirichlet.mle(Y[fall,])

## Get the corresponding alpha values, we round
## so we can report the values in the paper
## and this prevent any 0 probabilities (or close)
alpha.spring <- round(fit.sprg$alpha, 1)
alpha.summer <- round(fit.summ$alpha, 1)
alpha.fall <- round(fit.fall$alpha, 1)

noChangeDataFun <- function(n=63) {
  Y.sim <- matrix(nrow=n, ncol=4)
  Y.sim[sprg,] <- rdirichlet(n=n/3, alpha=alpha.spring)
  Y.sim[summ,] <- rdirichlet(n=n/3, alpha=alpha.summer)
  Y.sim[fall,] <- rdirichlet(n=n/3, alpha=alpha.fall)
  Y.sim
}

## A list of 200 datasets each of length 63 from the Dirichlet
set.seed(314159)
noChangeData <- lapply(rep(63, 200), noChangeDataFun)

## Now we generate data with a shift in location
##
## Here we base values on the first 7 years of our data
## and the last 7 years of data, accounting for seasonality

sprg1 <- seq(1,19,3)
summ1 <- seq(2,20,3)
fall1 <- seq(3,21,3)
sprg2 <- seq(43,61,3)
summ2 <- seq(44,62,3)
fall2 <- seq(45,63,3)

fit.sprg1 <- dirichlet.mle(Y[sprg1,])
fit.sprg2 <- dirichlet.mle(Y[sprg2,])
fit.summ1 <- dirichlet.mle(Y[summ1,])
fit.summ2 <- dirichlet.mle(Y[summ2,])
fit.fall1 <- dirichlet.mle(Y[fall1,])
fit.fall2 <- dirichlet.mle(Y[fall2,])

## The spread parameter, alpha0 or tau is based on the full data
## because we are not letting it change, fetch that from our
## fits on the full data
tau.spring <- round(fit.sprg$alpha0, 0)
tau.summer <- round(fit.summ$alpha0, 0)
tau.fall <- round(fit.fall$alpha0, 0)

## The theta values are from each the regimes (first or last 7 years)
theta.spring1 <- round(fit.sprg1$xsi,2)
theta.spring2 <- round(fit.sprg2$xsi,2)
theta.summer1 <- round(fit.summ1$xsi,2)
theta.summer2 <- round(fit.summ2$xsi,2)
theta.fall1 <- round(fit.fall1$xsi,2)
theta.fall2 <- round(fit.fall2$xsi,2)

changeLocationFun <- function(n=63, chpt=30) {
  Y.sim <- matrix(nrow=n, ncol=4, 0)
  sprg <- seq(1, 61, 3)
  summ <- seq(2, 62, 3)
  fall <- seq(3, 63, 3)
  Y.sim[sprg[(sprg<=chpt)],] <- rdirichlet(n=sum(sprg<=chpt), alpha=tau.spring*theta.spring1)
  Y.sim[sprg[(sprg>chpt)],] <- rdirichlet(n=sum(sprg>chpt), alpha=tau.spring*theta.spring2)
  Y.sim[summ[(summ<=chpt)],] <- rdirichlet(n=sum(summ<=chpt), alpha=tau.summer*theta.summer1)
  Y.sim[summ[(summ>chpt)],] <- rdirichlet(n=sum(summ>chpt), alpha=tau.summer*theta.summer2)
  Y.sim[fall[(fall<=chpt)],] <- rdirichlet(n=sum(fall<=chpt), alpha=tau.fall*theta.fall1)
  Y.sim[fall[(fall>chpt)],] <- rdirichlet(n=sum(fall>chpt), alpha=tau.fall*theta.fall2)
  Y.sim
}
## A list of 200 time series on simplex of length 63
set.seed(314159)
locationChangeData <- lapply(rep(63, 200), changeLocationFun)

############
## Similar story to the above but we let the
## spread parameter shift based on the regime (first or last 7 years),
## rather than the location which stays constant in time
## here we let the spread/scale change.  It is difficult to detect
## a shift in spread, so here we are letting it vary quite a bit
## An underlying scale shift of 12 units

tau.spring1.ex <- tau.spring + 4.5
tau.spring2.ex <- tau.spring - 4.5

tau.summer1.ex <- tau.summer - 4.5
tau.summer2.ex <- tau.summer + 4.5

tau.fall1.ex <- tau.fall - 4.5
tau.fall2.ex <- tau.fall + 4.5

tau.spring1 <- round(fit.sprg1$alpha0, 0)
tau.spring2 <- round(fit.sprg2$alpha0, 0)

tau.summer1 <- round(fit.summ1$alpha0, 0)
tau.summer2 <- round(fit.summ2$alpha0, 0)

tau.fall1 <- round(fit.fall1$alpha0,0)
tau.fall2 <- round(fit.fall2$alpha0,0)

theta.spring <- round(fit.sprg$xsi, 2)
theta.summer <- round(fit.summ$xsi, 2)
theta.fall <- round(fit.fall$xsi, 2)

###### Uncomment for the larger scale simulations
tau.spring1 <- tau.spring1.ex
tau.spring2 <- tau.spring2.ex
tau.summer1 <- tau.summer1.ex
tau.summer2 <- tau.summer2.ex
tau.fall1 <- tau.fall1.ex
tau.fall2 <- tau.fall2.ex

changeSpreadFun <- function(n=63, chpt=30) {
  Y.sim <- matrix(nrow=n, ncol=4, 0)
  sprg <- seq(1, 61, 3)
  summ <- seq(2, 62, 3)
  fall <- seq(3, 63, 3)
  Y.sim[sprg[(sprg<=chpt)],] <- rdirichlet(n=sum(sprg<=chpt), alpha=tau.spring1*theta.spring)
  Y.sim[sprg[(sprg>chpt)],] <- rdirichlet(n=sum(sprg>chpt), alpha=tau.spring2*theta.spring)
  Y.sim[summ[(summ<=chpt)],] <- rdirichlet(n=sum(summ<=chpt), alpha=tau.summer1*theta.summer)
  Y.sim[summ[(summ>chpt)],] <- rdirichlet(n=sum(summ>chpt), alpha=tau.summer2*theta.summer)
  Y.sim[fall[(fall<=chpt)],] <- rdirichlet(n=sum(fall<=chpt), alpha=tau.fall1*theta.fall)
  Y.sim[fall[(fall>chpt)],] <- rdirichlet(n=sum(fall>chpt), alpha=tau.fall2*theta.fall)
  Y.sim
}

set.seed(314159)
spreadChangeData <- lapply(rep(63, 200), changeSpreadFun)


##########################################
## Multiple change point problem
## At chpt1 there is a shift in location
## At chpt2 there is a shift in scale
changeLocSpreadMultCPFun <- function(n=63, chpt1=21, chpt2=42) {
  Y.sim <- matrix(nrow=n, ncol=4, 0)
  sprg <- seq(1, 61, 3)
  summ <- seq(2, 62, 3)
  fall <- seq(3, 63, 3)
  
  Y.sim[sprg[(sprg<=chpt1)],] <- rdirichlet(n=sum(sprg<=chpt1), alpha=tau.spring1*theta.spring1)
  Y.sim[sprg[(sprg>chpt1)&(sprg<=chpt2)],] <- rdirichlet(n=sum((sprg>chpt1)&(sprg<=chpt2)), alpha=tau.spring1*theta.spring2)
  Y.sim[sprg[(sprg>chpt2)],] <- rdirichlet(n=sum(sprg>chpt2), alpha=tau.spring2*theta.spring2)

  Y.sim[summ[(summ<=chpt1)],] <- rdirichlet(n=sum(summ<=chpt1), alpha=tau.summer1*theta.summer1)
  Y.sim[summ[(summ>chpt1)&(summ<=chpt2)],] <- rdirichlet(n=sum((summ>chpt1)&(summ<=chpt2)), alpha=tau.summer1*theta.summer2)
  Y.sim[summ[(summ>chpt2)],] <- rdirichlet(n=sum(summ>chpt2), alpha=tau.summer2*theta.summer2)

  Y.sim[fall[(fall<=chpt1)],] <- rdirichlet(n=sum(fall<=chpt1), alpha=tau.fall1*theta.fall1)
  Y.sim[fall[(fall>chpt1)&(fall<=chpt2)],] <- rdirichlet(n=sum((fall>chpt1)&(fall<=chpt2)), alpha=tau.fall1*theta.fall2)
  Y.sim[fall[(fall>chpt2)],] <- rdirichlet(n=sum(fall>chpt2), alpha=tau.fall2*theta.fall2)
  Y.sim
}

set.seed(314159)
locationSpreadChangeData <- lapply(rep(63, 200), changeLocSpreadMultCPFun)



##############################
## Data with covariate
plankton.data <- read.csv("biomass3.csv")
sprg <- seq(1,61,3)
summ <- seq(2,62,3)
fall <- seq(3,63,3)

summary(plankton.data[sprg,]$Epi_Temp)
sd(plankton.data[sprg,]$Epi_Temp)

summary(plankton.data[summ,]$Epi_Temp)
sd(plankton.data[summ,]$Epi_Temp)

summary(plankton.data[fall,]$Epi_Temp)
sd(plankton.data[fall,]$Epi_Temp)


covariateFun <- function(n=63) {
  Y.sim <- matrix(nrow=n, ncol=4, 0)
  X.sim <- matrix(nrow=n, ncol=1)
  sprg <- seq(1, 61, 3)
  summ <- seq(2, 62, 3)
  fall <- seq(3, 63, 3)
  
  X.sim[sprg,] <- rnorm(n=length(sprg), mean=16, sd=1.5)
  X.sim[summ,] <- rnorm(n=length(summ), mean=27, sd=1.5)
  X.sim[fall,] <- rnorm(n=length(fall), mean=21, sd=1.5)

  eta11 <- exp( -3   +  0.20*X.sim )
  eta12 <- exp(0.8  + -0.02*X.sim)
  eta13 <- exp(1.6  + -0.20*X.sim)
  eta14 <- exp(0.7 +  0.05*X.sim)
  
  eta21 <- exp( -1  + 0.20*X.sim)
  eta22 <- exp(0.8 + -0.05*X.sim)
  eta23 <- exp(1.4 + -0.05*X.sim)
  eta24 <- exp(0.3 + 0.10*X.sim)
  
  loc1 <- cbind(eta11, eta12, eta13, eta14)
  loc1.sum <- apply(loc1, 1, sum)
  loc1 <- apply(loc1, 2, function(x) { x/loc1.sum } )
  
  loc2 <- cbind(eta21, eta22, eta23, eta24)
  loc2.sum <- apply(loc2, 1, sum)
  loc2 <- apply(loc2, 2, function(x) {x/loc2.sum} )
  
  scal1 <- exp(-0.2 + 0.10*X.sim)
  scal2 <- exp(-0.15 + 0.15*X.sim)

  shape1 <- apply(loc1, 2, function(x) {x*scal1 } )
  shape2 <- apply(loc2, 2, function(x) {x*scal2 } )
  shape <- rbind(shape1[1:30,], shape2[31:63,])
  
  Y.sim <- rdirichlet(n=n, alpha=shape)
  cbind(Y.sim, X.sim)
}

set.seed(314159)
covariateData <- lapply(rep(63, 200), covariateFun)



##################################
## Data with a single change in location
## and scale, but done in such a way the 
## shape does not change too much

changeLocationSpreadFun <- function(n=63, chpt=30) {
  Y.sim <- matrix(nrow=n, ncol=4, 0)

  loc1 <- c(0.3, 0.3, 0.2, 0.2)
  scal1 <- 6
  loc1*scal1
  
  loc2 <- c(0.25, 0.25, 0.25, 0.25)
  scal2 <- 4
  loc2*scal2
  
  Y.sim[1:31,] <- rdirichlet(31, alpha=loc1*scal1)
  Y.sim[32:n,] <- rdirichlet(32, alpha=loc2*scal2)
  Y.sim
}

set.seed(314159)
locationSpreadChangeData <- lapply(rep(63, 200), changeLocationSpreadFun)


###
## Now save all these data for future use
##
save(uniformNoChangeData, file="SimulatedUniformNoChange.RData")
save(noChangeData, file="SimulatedNoChange.RData")
save(locationChangeData, file="SimulatedLocationChange.RData")
save(spreadChangeData, file="SimulatedSpreadChange.RData")

save(covariateData, file="SimulatedCovariate.RData")
save(locationSpreadChangeData, file="SimulatedLocationSpreadChange.RData")

## We will also save all the parameters so we can report them in the paper
simulatedDataParameters <- list(alpha.spring=alpha.spring, 
              alpha.summer=alpha.summer, 
              alpha.spring=alpha.spring,
              theta.spring = theta.spring,
              theta.summer = theta.summer,
              theta.fall = theta.fall,
              tau.spring = tau.spring,
              tau.summer = tau.summer,
              tau.fall = tau.fall,
              
              theta.spring1 = theta.spring1,
              theta.summer1 = theta.summer1,
              theta.fall1 = theta.fall1,
              tau.spring1 = tau.spring1,
              tau.summer1 = tau.summer1,
              tau.fall1 = tau.fall1,
              tau.spring1.ex = tau.spring1.ex,
              tau.summer1.ex = tau.summer1.ex,
              tau.fall1.ex = tau.fall1.ex,
              
              theta.spring2 = theta.spring2,
              theta.summer2 = theta.summer2,
              theta.fall2 = theta.fall2,
              tau.spring2 = tau.spring2,
              tau.summer2 = tau.summer2,
              tau.fall2 = tau.fall2,
              tau.spring2.ex = tau.spring2.ex,
              tau.summer2.ex = tau.summer2.ex,
              tau.fall2.ex = tau.fall2.ex
              
              )

save(simulatedDataParameters, file="SimulatedDataParameters.RData")

