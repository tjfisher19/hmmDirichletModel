######################
## Code that runs the simulated data
##  on the nonparameteric test

source("nonparametricChangePoint.R")

#################################
## First chunk runs the nonparametric test on the 
## simulated datasets.

load("convergenceTable.RData")

load("SimulatedNoChange.RData")

cpDistNoChange <- sapply(ind.picked$NoChangeLocSpreadModel, function(x) {getNonParametricChangePoint(noChangeData[[x]], b=1) } )
ecpNoChange <- lapply(ind.picked$NoChangeLocSpreadModel, function(x) { getEcpTestChangePoint(noChangeData[[x]])})

mean(cpDistNoChange[1,]<0.05)
mean(unlist(lapply(ecpNoChange, function(x) {length(x)>2})))

rm(noChangeData)




load("SimulatedLocationChange.RData")

cpDistLocationChange <- sapply(ind.picked$LocChangeLocSpreadModel, function(x) {getNonParametricChangePoint(locationChangeData[[x]], b=1) } )
ecpLocationChange <- lapply(ind.picked$LocChangeLocSpreadModel, function(x) { getEcpTestChangePoint(locationChangeData[[x]])})

mean(cpDistLocationChange[1,]<0.05)
mean(unlist(lapply(ecpLocationChange, function(x) {length(x)>2})))




load("SimulatedSpreadChange.RData")

cpDistSpreadExChange <- sapply(ind.picked$SpreadExLocSpreadModel, function(x) {getNonParametricChangePoint(spreadChangeData[[x]], b=1) } )
ecpSpreadExChange <- lapply(ind.picked$SpreadExLocSpreadModel, function(x) { getEcpTestChangePoint(spreadChangeData[[x]])})

mean(cpDistSpreadExChange[1,]<0.05)
mean(unlist(lapply(ecpSpreadExChange, function(x) {length(x)>2})))


load("tom.RData")
load("Simulated3stateLocationSpreadChange.RData")

cpDistLocSpread3state <- sapply(ind.picked$locScale3State, function(x) {getNonParametricMultChangePoint(locationSpreadChangeData3[[x]], b=1) } )
ecpLocSpread3State <- lapply(ind.picked$locScale3State, function(x) {getEcpTestChangePoint(locationSpreadChangeData3[[x]]) } )

mean(sapply(1:100, function(x) {length(cpDistLocSpread3state[[x]]) } )==1)
mean(sapply(1:100, function(x) {length(cpDistLocSpread3state[[x]]) } )==2)

mean(sapply(1:100, function(x) {length(ecpLocSpread3State[[x]]) } )==3)
mean(sapply(1:100, function(x) {length(ecpLocSpread3State[[x]]) } )==4)


load("convergenceHardDetect.RData")
load("SimulatedLocationSpreadChange.RData")

cpDistLocationSpreadHardChange <- sapply(ind.picked$hardDetect, function(x) {getNonParametricChangePoint(locationSpreadChangeData[[x]], b=1) } )
ecpLocationSpreadHardChange <- lapply(ind.picked$hardDetect, function(x) { getEcpTestChangePoint(locationSpreadChangeData[[x]])})

mean(cpDistLocationSpreadHardChange[1,]<0.05)
mean(unlist(lapply(ecpLocationSpreadHardChange, function(x) {length(x)>2})))


save(cpDistNoChange, 
     ecpNoChange,
     cpDistLocationChange, 
     ecpLocationChange,
     cpDistSpreadExChange, 
     ecpSpreadExChange,
     cpDistLocSpread3state,
     ecpLocSpread3State,
     cpDistLocationSpreadHardChange,
     ecpLocationSpreadHardChange,
     file="nonparametricTestResults.RData")

