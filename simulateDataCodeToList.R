##########################################
##
## This converts the original simulated data
## into a list of list so they can
## easily be incorporated into RStan
##
## This was done to improve rstan efficiency
## on the OSC cluster.
##
## This code basically takes the simulation data
## and adds other parameters necessary for the rstan code
## It puts it all into a list() so the parLapply() functions
## on the cluster code well.

load("SimulatedUniformNoChange.RData")

convertToList <- function(ind) {
    within(list(), {
      T <- 63
      M <- 2
      K <- 4
      P <- 1
      y <- uniformNoChangeData[[ind]]
      X <- matrix(nrow=dim(uniformNoChangeData[[ind]])[1], ncol=1, 1)
      rho_a <- 9.5
      rho_b <- 0.5
    })
}

simulatedData <- lapply(1:length(uniformNoChangeData), convertToList)

save(simulatedData, file="ListedUniformNoChange.RData")




load("SimulatedNoChange.RData")

convertToList <- function(ind) {
  within(list(), {
    T <- 63
    M <- 2
    K <- 4
    P <- 3
    y <- noChangeData[[ind]]
    X <- as.matrix(model.matrix(~as.factor(rep(1:3, 21))) )
    rho_a <- 9.5
    rho_b <- 0.5
  })
}

simulatedData <- lapply(1:length(noChangeData), convertToList)

save(simulatedData, file="ListedNoChange.RData")




load("SimulatedLocationChange.RData")

convertToList <- function(ind) {
  within(list(), {
    T <- 63
    M <- 2
    K <- 4
    P <- 3
    y <- locationChangeData[[ind]]
    X <- as.matrix(model.matrix(~as.factor(rep(1:3, 21))) )
    rho_a <- 9.5
    rho_b <- 0.5
  })
}

simulatedData <- lapply(1:length(locationChangeData), convertToList)

save(simulatedData, file="ListedLocationChange.RData")



load("SimulatedSpreadChange.RData")

convertToList <- function(ind) {
  within(list(), {
    T <- 63
    M <- 2
    K <- 4
    P <- 3
    y <- spreadChangeData[[ind]]
    X <- as.matrix(model.matrix(~as.factor(rep(1:3, 21))) )
    rho_a <- 9.5
    rho_b <- 0.5
    b_loc_sig <- 2.0
    b_scale_sig <- 2.0
  })
}

simulatedData <- lapply(1:length(spreadChangeData), convertToList)

save(simulatedData, file="ListedSpreadChange.RData")





load("SimulatedCovariate.RData")

convertToList <- function(ind) {
  within(list(), {
    T <- 63
    M <- 2
    K <- 4
    P <- 2
    y <- covariateData[[ind]][,1:4]
    X <- as.matrix(model.matrix(~covariateData[[ind]][,5]) )
    rho_a <- 9.5
    rho_b <- 0.5
    b_loc_sig <- 2.0
    b_scale_sig <- 2.0
  })
}

simulatedData <- lapply(1:length(covariateData), convertToList)

save(simulatedData, file="ListedCovariateChange.RData")

convertToList <- function(ind) {
  within(list(), {
    T <- 63
    M <- 2
    K <- 4
    P <- 3
    y <- covariateData[[ind]][,1:4]
    X <- as.matrix(model.matrix(~as.factor(rep(1:3, T/3))) )
    rho_a <- 9.5
    rho_b <- 0.5
    b_loc_sig <- 2.0
    b_scale_sig <- 2.0
  })
}

simulatedData <- lapply(1:length(covariateData), convertToList)

save(simulatedData, file="ListedCovariateChangeSeasons.RData")


load("SimulatedLocationSpreadChange.RData")

convertToList <- function(ind) {
  within(list(), {
    T <- 63
    M <- 2
    K <- 4
    P <- 1
    y <- locationSpreadChangeData[[ind]][,1:4]
    X <- as.matrix(rep(1, T) )
    rho_a <- 8.0
    rho_b <- 2.0
    b_loc_sig <- 2.0
    b_scale_sig <- 2.0
  })
}

simulatedData <- lapply(1:length(locationSpreadChangeData), convertToList)

save(simulatedData, file="ListedLocationSpreadChange.RData")

