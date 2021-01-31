##############################
## Real phytoplankton data
##  Here we run the nonparametric &
##  Dirichlet likelihood method.

source("dirichletPermutationTest.R")
source("nonparametricChangePoint.R")

plankton.data <- read.csv("biomass3.csv")
Y <- plankton.data %>%
  dplyr::select(Prop.BG, Prop.Green, Prop.Flag, Prop.Diatom)

transData <- as.data.frame(log(Y[,1:3]/Y[,4]))
transData$Season <- rep(1:3, 21)
transData.summary <- transData %>%
  group_by(Season) %>%
  summarize(Mean.V1 = mean(Prop.BG),
            Mean.V2 = mean(Prop.Green),
            Mean.V3 = mean(Prop.Flag))

transData <- left_join(transData, transData.summary, by="Season") %>%
  mutate(V1 = Prop.BG - Mean.V1,
         V2 = Prop.Green - Mean.V2,
         V3 = Prop.Flag - Mean.V3) %>%
  dplyr::select(V1, V2, V3)

npOut1 <- cpDist(as.matrix(transData), b=1)
names(Y) <- c("V1", "V2", "V3", "V4")
npOut2 <- getNonParametricMultChangePoint(Y)
set.seed(2020)
ecpOut1 <- e.divisive(as.matrix(transData), min.size=18)
set.seed(2020)
dirPermTest <- parametricDirichletMultCP(as.matrix(Y) )

