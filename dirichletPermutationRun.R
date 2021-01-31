###################
## Code that calls the Dirichlet likelihood
##  permutation test on the simulated data


library(parallel)

options(mc.cores = 8, width=150)

source("dirichletPermutationTest.R")

cl <- makeForkCluster()

# load("SimulatedNoChange.RData")
# load("summary/convergenceTable.RData")
# 
# cpDirNoChange <- parLapply(cl, noChangeData[ind.picked$NoChangeLocSpreadModel], parametricDirichletTest)
# # cpDirNoChange <- bind_rows(cpDirNoChange)
# 
# rm(noChangeData)
# 
# load("SimulatedLocationChange.RData")
# 
# cpDirLocChange <- parLapply(cl, locationChangeData[ind.picked$LocChangeLocSpreadModel], parametricDirichletTest)
# # cpDirLocChange <- bind_rows(cpDirNoChange)
# 
# rm(locationChangeData)
# 
# load("SimulatedSpreadChangeBackup.RData")
# 
# 
# load("SimulatedSpreadChange.RData")
# 
# cpDirSpreadExChange <- parLapply(cl, spreadChangeData[ind.picked$SpreadChangeLocSpreadModel], parametricDirichletTest)
# # cpDirSpreadExChange <- bind_rows(cpDirSpreadExChange)
# 
# rm(spreadChangeData)
# 
# save(cpDirNoChange, 
#      cpDirLocChange, 
#      cpDirSpreadExChange,
#      file="dirichletPermutedTestResults.RData")



load("tom.RData")
load("SimulatedLocationSpreadChange.RData")

cpDirLocationSpreadChange <- parLapply(cl, locationSpreadChangeData[ind.picked$locScale3State], parametricDirichletMultCP)


stopCluster(cl)

save(cpDirLocationSpreadChange,
     file="dirichletPermutedTestMultCPResults.RData")




