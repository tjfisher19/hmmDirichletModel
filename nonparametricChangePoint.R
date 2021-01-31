########################
## Our implementation of the 
## Multivariate nonparametric tests for a change in distribution
## and the Bayesian segmentation test
##
## For all these tests we need to transform the data
## 
## We take the 4-dimensional compositional data and perform
## a log ratio conditioning on the 4th element, thus transforming
## from the simplex of 4 dimensions into the trivariate real
## number system. We then run the test with a specified bandwidth b
##
##
## For the npcp test
##
## The bandwidth b determines the modeling of serial dependence
##  A value of b=1 assumes independence. Given the time series nature
##    of our data, this is likely incorrect, but it also results in the
##    highest levels of power for this test, so we report it in the paper
##  If you let b=NULL, then the algorithm will estimate an optimal value
##    of b using the bOptEmpProc() function. This is a lot slower
##   
##  In the function, getNonParametricMultChangePoint(), we perform the 
##    algorithm outlined in section 2.4 in Prabuchandran et al. (2019)
##    but apply the nonparametric test. 
##
## For the energy divisive (ecp)
##
##  Also below is code that implements the Energy Divisive test.
##   That test can work for a general distribution but assumes
##   the data is in R^d, not a simplex. We perform the same
##   transformation from 4 dimensions into 3.
##
##   
library(npcp)
library(ecp)
library(dplyr)

getNonParametricChangePoint <- function(data, b=1) {
  transData <- as.data.frame(log(data[,1:3]/data[,4]))
  transData$Season <- rep(1:3, 21)
  transData.summary <- transData %>%
    group_by(Season) %>%
    summarize(Mean.V1 = mean(V1),
              Mean.V2 = mean(V2),
              Mean.V3 = mean(V3))
  
  transData <- left_join(transData, transData.summary, by="Season") %>%
    mutate(V1 = V1 - Mean.V1,
           V2 = V2 - Mean.V2,
           V3 = V3 - Mean.V3) %>%
    dplyr::select(V1, V2, V3)
  
  tmp <- cpDist(as.matrix(transData), b=b)
  c(tmp$p.value, which.max(tmp$cvm))
}


getNonParametricMultChangePoint <- function(data, b=1) {
  transData <- as.data.frame(log(data[,1:3]/data[,4]))
  transData$Season <- rep(1:3, 21)
  transData.summary <- transData %>%
    group_by(Season) %>%
    summarize(Mean.V1 = mean(V1),
              Mean.V2 = mean(V2),
              Mean.V3 = mean(V3))
  
  transData <- left_join(transData, transData.summary, by="Season") %>%
    mutate(V1 = V1 - Mean.V1,
           V2 = V2 - Mean.V2,
           V3 = V3 - Mean.V3) %>%
    dplyr::select(V1, V2, V3)
  
  ## Using the algorithm in Prabuchandran et al. (2019)
  ##  Start off with L terms (L<n), look for changepoint
  ## If you find one at time point c, then look again
  ## with a series from c+1, c+2, ..., C+L
  ## For convenience we will use L = 30
  ## This should help the algorithm detect the multiple
  ## change points
  n <- dim(transData)[1]
  L0 <- 30
  L <- L0
  batch <- 10
  begin <- 1
  cpt <- NULL
  run <- TRUE
  while(run) {
    if(L >= n) {
      L <- n
      run=FALSE
    }
    tmp <- cpDist(as.matrix(transData)[begin:L,], b=b)
    if(tmp$p.value<0.05) {
      cpt <- c(cpt, which.max(tmp$cvm))
      begin <- L+1
      L <- L+L0

    } else {
      L <- L + batch
    }
  }
  cpt
}

####
## Since our real data and simulated are only based on 63 observations
## with 3 measuremetns per year. We set the minimum segmentation length
## to be 18 (6 years).  
getEcpTestChangePoint <- function(data) {
  transData <- as.data.frame(log(data[,1:3]/data[,4]))
  transData$Season <- rep(1:3, 21)
  transData.summary <- transData %>%
    group_by(Season) %>%
    summarize(Mean.V1 = mean(V1),
              Mean.V2 = mean(V2),
              Mean.V3 = mean(V3))
  
  transData <- left_join(transData, transData.summary, by="Season") %>%
    mutate(V1 = V1 - Mean.V1,
           V2 = V2 - Mean.V2,
           V3 = V3 - Mean.V3) %>%
    dplyr::select(V1, V2, V3)
  
  e.divisive(as.matrix(transData), min.size=18)$estimates
}

