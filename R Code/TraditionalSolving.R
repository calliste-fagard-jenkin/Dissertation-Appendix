source("SimulateData.R")
source("piecemealFunctions.R")
source("PlottingFunctions.R")
library(plot.matrix) # for heatmaps of transition

# Generate some data:
set.seed(102938)
tsteps <- 120
SimmedData <- simulateIBM(n=50, t=tsteps,
                          # set survival details:
                          survFunc = linLogit, survPars = c(2.26, 0.23),
                          # set growth details:
                          growthSamp = sampleDTN,
                          growthPars = c(0, 1, log(0.1), 0, 10),
                          # set reproduction details:
                          reprFunc = linLogit, reprPars = c(-3.5, 0.25),
                          # set offspring number and size distribution details:
                          offNumSamp = sampleOffNum, offNumPars = log(4),
                          offSizeSamp = sampleDTN,
                          offSizePars = c(1, 0.6, log(0.1), 0, 10),
                          Schild=0.6,
                          # set other miscelaneous parameters:
                          Start=3, thresh=1000, OneGend = TRUE, popPrint=TRUE)
                       
# remove the observations we can't use for survival fitting:
solveDF <- subset(SimmedData, !is.na(prev.size))

# Make a DF with only individuals that had the chance to reproduce:
solveDF2 <-  subset(SimmedData, !is.na(SimmedData$reproduced))

# Make a DF with only the individuals when they were born:
solveDF3 <- subset(SimmedData, !is.na(SimmedData$parent.size))

# Make a DF to estimate offspring per parent:
solveDF4 <- subset(SimmedData, !is.na(SimmedData$off.born))

# Make a DF to estimate offspring survival:
deadChildren <- sum(solveDF4$off.born) - sum(solveDF4$off.survived)
solveDF5 <- data.frame(survived=c(solveDF3$survived, rep(0, deadChildren)))

############## PIECEMEAL ESTIMATION OF DEMOGRAPHIC PARAMETERS ##################
survivalSol <- glm(survived ~ prev.size, family = binomial, data = solveDF)
growthSol <- lm(size ~ prev.size, data = solveDF)
reproductionSol <- glm(reproduced ~ size, family = binomial, data=solveDF2)
offspringSizeSol <- lm(size ~ parent.size, data = solveDF3)
offCountSol <- glm(off.born ~ 1, family = poisson, data = solveDF4)
offspringSurvSol <- glm(survived ~ 1, family=binomial, data=solveDF5)

# survivalSol %>% coef %>% print
# growthSol %>% coef %>% print
# reproductionSol %>% coef %>% print
# offspringSizeSol %>% coef %>% print
# offCountSol %>% coef %>% exp %>% print
# offspringSurvSol %>% coef %>% plogis %>% print

######################### IMPLIMENTING THE IPM #################################

# Calculate the projection matrix using estimated parameters:
A <- kernelOneVar(m=100,
                  # growth:
                  growthFunc = doublyTruncatedNormal,
                  growthPars = c(coef(growthSol),
                                 growthSol %>% summary %>% `$`(sigma) %>% log,
                                 0, 10),
                  # survival:
                  survFunc = linLogit,
                  survPars = coef(survivalSol),
                  # reproduction:
                  repFunc = linLogit,
                  repPars = coef(reproductionSol),
                  # offspring number per birth:
                  offNum = offCountSol %>% coef %>% exp,
                  # offpspring sizes:
                  offSizeFunc = doublyTruncatedNormal,
                  offSizePars = c(coef(offspringSizeSol), offspringSizeSol %>%
                                    summary %>% `$`(sigma) %>% log, 0, 10),
                  # child survival rate:
                  childSurv = offspringSurvSol %>% coef %>% plogis,
                  # Other parameters:
                  halfPop = TRUE, L=0, U=10)


plot(A)
A %>% eigen %>% `$`(values) %>% `[`(1) %>% Re
