source('AdaptivePMCMC.R')
source('simulateData.R')
source('ModelSpecIPM.R')
source('piecemealFunctions.R')
source('SMC.R')
library(dissPackage3)
library(xtable)
library(mcmcse)

########################### SPECIFYING THE PRIORS ##############################

# Setup the prior names:
priorsIPM <- rep("dunif", 14)

# Setup the prior hyper par values:
flatPriors <- list(
  # Survival function intercept:
  sFuncInt=list(min=-50, max=50),
  # Survival function gradient:
  sFuncGra=list(min=-50, max=50),
  # Growth function intercept:
  gFuncInt=list(min=-50, max=50),
  # Growth function gradient:
  gFuncGra=list(min=-50, max=50),
  # Growth function sd:
  gFuncSD=list(min=-50, max=50),
  # Reproduction function intercept:
  rFuncInt=list(min=-50, max=50),
  # Reproduction function gradient:
  rFuncInt=list(min=-50, max=50),
  # Number of offpsring per litter;
  offNum=list(min=-50, max=50),
  # Childsize intercept:
  osFuncInt=list(min=-50, max=50),
  # Childsize gradient:
  osFuncGra=list(min=-50, max=50),
  # Childsize sd:
  osFuncSD=list(min=-50, max=50),
  # Child survival:
  Schild=list(min=-50, max=50),
  # The probability of detection of an animal:
  obsProb=list(min=-50, max=50)
)

skeleton = list(
  survPars = rep(NA, 2),
  growthPars = rep(NA, 3),
  reprPars = rep(NA, 2),
  offNumPars = NA,
  offSizePars = rep(NA, 3),
  Schild = NA,
  obsProbPar = NA
)

########## CALCULATE SOME VALUES BEFORE WE CAN MAKE THE IPMLTP OBJECT ##########

simPars <- list(n=100, t=5000,
                # set survival details:
                survFunc = linLogit, survPars = c(-9.65, 3.77),
                # set growth details:
                growthSamp = sampleDTN,
                growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
                # set reproduction details:
                reprFunc = linLogit, reprPars = c(-7.23, 2.6),
                # set offspring number and size distribution details:
                offNum=1, offSizeSamp = sampleDTN,
                offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
                # Child survival probability:
                Schild=0.873,
                # set other miscelaneous parameters:
                Start = 2.7, thresh=700, OneGend = TRUE,
                popPrint = T, verbose=F)

# Define the breaks for the 4 group specification:
breaks <- seq(1.5, 3.55, l=6)[-2]

# Get the prior initial size distribution:
priorLoopN <- 100
priorCounts <- rep(0, length(breaks)-1)
priorStartSize <- rep(NA, priorLoopN)

set.seed(201020)
for (i in 1:priorLoopN){
  # Generate the data:
  priorData <- do.call(simulateIBM, simPars)  
  max.prior <- priorData$census.number %>% max
  # Extract the counts in each size class for the first census:
  loop.counts <- subset(priorData, priorData$census.number==(max.prior-5)) %>%
    `$`(size) %>% na.omit %>% vectorToCounts(breaks = breaks)
  # Update the previous vector:
  priorCounts %<>% `+`(loop.counts)
  priorStartSize[i] <- sum(loop.counts)
}

# Get the parameters for the initial size distribution:
priorProbs <- countsToProbs(priorCounts)
priorStartSize <- (sum(priorCounts)/priorLoopN) %>% ceiling
priorStartSizeSD <- sd(priorStartSize)
############################ GENERATE THE DATA #################################

# Get the observations:
set.seed(102010)
simmedData <- do.call(simulateIBM, simPars)
max.cens <- simmedData$census.number %>% max
Y <- getSizeDistns(simmedData, breaks)[,(max.cens-5):max.cens]

#################### CREATE THE LOGTARGETPARAMETERS OBJECT #####################

IPMLTP <- list(
  priorFunc = 'evalPriors',
  priors = priorsIPM,
  priorPars = flatPriors,
  skeleton = skeleton,
  survFunc = 'linLogit',
  growthSamp = 'sampleNorm',
  reprFunc = 'linLogit',
  offNumSamp = 'returnConstant',
  offSizeSamp = 'sampleNorm',
  oneSex = TRUE,
  mu = 'multnomMu',
  muPar = c(sum(Y[,1]), list(priorProbs)),
  b = 250, 
  Y = Y,
  obsProb = 'detectionNumObs',
  shift = qlogis(0.49),
  breaks = breaks
)

# saveRDS(IPMLTP, "IPMLTP")

########################### MAKE THE START VALUES ##############################

startValues <- list(
  survPars = c(-9.65, 3.77),
  growthPars = c(1.41, 0.56, log(0.08)),
  reprPars = c(-7.23, 2.6),
  offNumPars = 1,
  offSizePars = c(0.36, 0.71, log(0.16)),
  Schild = qlogis(0.873), 
  obsProbPar = 10 # not too close to 50 since this will hurt the chain
)

startValues %<>% unlist
SVjitter <- startValues + rnorm(length(startValues), 0, 0.5)

########################### SIMULATED ANNEALING ################################


set.seed(101010)
# get the hessian:
SANN <-  optim(startValues, logTargetIPM, logTargetPars = IPMLTP, returnNeg = T,
               method = "SANN", control = list(maxit = 10), printProp = T,
               hessian = T)
# saveRDS(SANN, "SANN")
# get the covariance of the proposal:
propCov <- solve(SANN$hessian)

set.seed(101010)
IPMLTP2 <- IPMLTP
IPMLTP2$b <- 1000
# get the MLEs by using more iterations:
MLEs <- optim(startValues, logTargetIPM, logTargetPars = IPMLTP2, returnNeg = T,
               method = "SANN", control = list(maxit = 1000), printProp = T)
# saveRDS(MLEs, "MLEs")
# Check some log lik values
logTargetIPM(startValues, logTargetPars = IPMLTP, returnNeg = T, printProp = F)
logTargetIPM(MLEs$par, logTargetPars = IPMLTP, returnNeg = T, printProp = F)
logTargetIPM(SVjitter, logTargetPars = IPMLTP, returnNeg = T, printProp = F)


######################### RUN THE MCMC CHAINS IN PARALLEL ######################

# high resolution plots:
plotNames <- c("log(pi)", "Si", "Sg", "Gi", "Gg", "log(Gs)", "Ri", "Rg", "ON",
               "OSi", "OSg", "log(OSs)", "qlogis(Schild)", "qlogis(ObsP)")

# Objects which the cluster needs to know about to run the sampling:
clNeeds = c('logTargetIPM', '%<>%', '%>%', 'sampleNorm', 'returnConstant',
            'detectionNumObs', 'evalPriors', 'sampleStateIPM', 'particleFilter',
            'multnomMu', 'linLogit', 'vectorisedSamplerIPM', 'rmvnorm') 

# An example of how chains are run, using a stored proposal covariance structure
# and then stored for iteration:
ungulateChainMult <- pMH(proposal = multvarNormProp, uFunc = NULL,
                         propPars = diag(length(startValues))/100,
                         #propPars = readRDS("ObsNumSigma4.1"),
                         #propPars = propCov,
                         lTarg = logTargetIPM, lTargPars = IPMLTP,
                         x0 = MLEs$par, itermax=1000, prntPars = T,
                         clNeeds = clNeeds, packages = "dissPackage3")

# plot(ungulateChainMult, cols=2:14, width=8.3*3, height=11.7*3,
#      cex=2, names=plotNames, filePath="")
# 
# saveRDS(cov(ungulateChainMult[[1]][10001:20000,-1]), "ObsNumSigma4.1")
# saveRDS(ungulateChainMult, 'SANNchain2')

####################### DIAGNOSTICS ON THE CHOSEN CHAIN ########################

# load the chain and get diagnostics:
chosenChain <- readRDS('bluewhale33')
getAcceptance(chosenChain)
grDiagnostic(chosenChain)

# make plots for the first three parameters as figures for the writeup:
chainFirst3 <- rep(NA, 3) %>% list
for (j in 1:4) chainFirst3[[j]] <- chosenChain[[j]][,-(4:14)] 
class(chainFirst3) <- 'pMCMCoutput'
F3plotNames <- plotNames[1:3]
grDiagnostic(chainFirst3)
plot(chainFirst3, cols = 2:3, width=8.3*3*(0.3), height=11.7*1.5*(0.3),
     cex=1, names=F3plotNames, filePath="chainFirst3_")

# plot chain before thinning:
plot(chosenChain, cols=2:14, width=8.3*3, height=11.7*3,
     cex=2, names=plotNames, filePath="")

# plot chain after thinning:
thinnedChosen <- thinMCMC(chosenChain, alpha = 0.1)
plot(thinnedChosen, cols=2:14, width=8.3*3, height=11.7*3,
     cex=2, names=plotNames, filePath="")

# sub-plot for the write-up:
plot(thinnedChosen, cols=c(2,4,7,10,13,14), width=8, height=8,
     cex=2, names=plotNames[c(1,2,4,7,10,13,14)], filePath="")

# effective sample sizes for the chains:
multiESS(thinnedChosen[[1]][,-1]) ; multiESS(thinnedChosen[[2]][,-1]) 
multiESS(thinnedChosen[[3]][,-1]) ; multiESS(thinnedChosen[[4]][,-1])

# combine the thinned chain samples together:
combinedChosen <- thinnedChosen[[1]]
for (i in 2:4) combinedChosen <- rbind(combinedChosen, thinnedChosen[[i]])

# apply link functions:
returnSelf <- function(x) x
links <- rep('returnSelf', 13)
links[c(5,11)] <- 'exp'
links[12:13] <- 'plogis'
for (i in 2:14) combinedChosen[,i] <- match.fun(links[i-1])(combinedChosen[,i])

# get the MAP parameters:
ind <- combinedChosen[,1] %>% which.max
MAP <-  combinedChosen[ind, -1]

# get the median, mean and 95% CIs:
means <- apply(combinedChosen[,-1], 2, mean)
medis <- apply(combinedChosen[,-1], 2, median)
lower <- apply(combinedChosen[,-1], 2, quantile, probs = 0.025)
upper <- apply(combinedChosen[,-1], 2, quantile, probs = 0.975)

# make a data.frame, and then tables for the write-up:
simulated <- c(-9.65, 3.77, 1.41, 0.56, 0.08, -7.23,
               2.6, 1, 0.36, 0.71, 0.16, 0.873, 1)

summaries <- data.frame(simulated = simulated, MAP = MAP, mean = means,
                        median = medis, lower = lower, upper = upper)

rownames(summaries) <- c("survival.i", "survival.g", "growth.i", "growth.g",
                         "growth.sd", "repr.i", "repr.g", "litter.size",
                         "off.size.i", "off.size.g", "off.size.sd",
                         "child.survival.p", "observation.p")

xtable(t(summaries)[,c(1:5, 13)])
xtable(t(summaries)[,-c(1:5, 13)])

# posterior growth rate distn from the samples:
cluster <- makeCluster(4)
CI <- posteriorGrowthRate(thinnedChosen,
                          IPMLTP, growthFunc = "normal", printRate = 10,
                          offSizeFunc = "normal", L = 1.5, U = 3.55, m = 500,
                          cluster = cluster, clNeeds = clNeeds)
