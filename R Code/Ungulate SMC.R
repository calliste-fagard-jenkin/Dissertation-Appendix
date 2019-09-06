source("SMC.R")
source("AdaptivePMCMC.R")
source("piecemealFunctions.R")
source("ModelSpecIPM.R")
source("simulateData.R")
source("plottingFunctions.R")

############################# SIMULATE THE DATA SET ############################

set.seed(102030)
# We simulate data from a constant size distribution. Projecting this forwards
# allows us to get a bette restimate of the starting distribution to simulate
# from:
sampleN=50
simmedData <- simulateIBM(n=sampleN, t=1000,
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
                          Start=3, thresh=5000, OneGend = TRUE, popPrint = T)

# Sample sampleN values from the final time distribution as an empirical
# estimate of the stable size distribution:
maxTime <- simmedData$census.number %>% max
sampleStart <- simmedData %>% subset(simmedData$census.number==maxTime) %>%
  `$`(size) %>% na.omit %>% sample(sampleN)

simmedData <- simulateIBM(n=sampleN, t=20,
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
                          Start = sampleStart, thresh=500, OneGend = TRUE,
                          popPrint = T, verbose=F)

# Uncomment the below code to keep only the last 5 censuses:
simmedData <- with(simmedData, # Keep only the last 5 censuses for stability
                   subset(simmedData, census.number>=(max(census.number)-5)))

# Reajust the census.number labels:
simmedData$census.number <- simmedData$census.number -
  min(simmedData$census.number) + 1

######################### MAKE SOME ILLUSTRATIVE PLOTS  ########################

breaks <- seq(1.5, 3.55, l=6)
sizeObservations <- getSizeDistns(simmedData, breaks=breaks)

# Produce a figure:
plot(log(sizeObservations), ylab="log count", xlab='Census number', main='')

# As a placeholder for the prior, use the distribution of observed sizes
# across all censuses:
muPs <- simmedData$size %>% na.omit %>% vectorToCounts(breaks) %>% countsToProbs


# Make a list of the parameters for the particle filter:
ungulatePars <- list(
  survFunc = 'linLogit', survPars = c(-9.65, 3.77),
  growthSamp = 'sampleDTN', growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
  reprFunc = 'linLogit', reprPars = c(-7.23, 2.6), 
  offSizeSamp = 'sampleDTN', offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
  offNumSamp = 'returnConstant', offNumPars = 1,
  Schild=0.873, oneSex = TRUE, breaks = breaks, shift=qlogis(0.44))

# Simulate a single time step:
do.call(sampleStateIPM, c(list(multnomMu(sampleN, muPs)), ungulatePars, verbose=T))

# Simulate multiple time steps:
obsPopSizes <- getSizeDistns(simmedData, breaks) %>% apply(2, sum) %>% print

plot(obsPopSizes, type='l', col=4, lwd=2, ylab='Population size',
     xlab='census number')

########################## USING THE PARTICLE FILTER ###########################

# Make a fake true population size vector:
popSizes <- sizeObservations %>% apply(2, sum) %>% `*`(1.2)
trueSizes <- makeFakeTrueSizeDist(muPs, popSizes)

# To illustrate the probability of our observations given the true pop size:
exampleObsProb <- multinomialObs(sizeObservations, trueSizes, 0.8)

# To illstrate underflow, attempt to normalise:
exampleObsProb/sum(exampleObsProb)
exampleObsProb - max(exampleObsProb)

trueStart <- sizeObservations[,1] %>% sum %>% `*`(1.2)

# Estimate the log likelihood at the known parameter values:
particleFilter(Y = sizeObservations,                             # Data
               mu = multnomMu, muPar = c(trueStart, list(muPs)), # Initial distn 
               vectorisedSamplerIPM, ungulatePars,               # SS info
               # Observation model specification:
               obsProb = detectionNumObs, obsProbPar = list(0.8),
               b = 100, returnW = F) # Particle filter specs

# about 4.8 to 5.2
#################### DOING AN ANALYSIS WITH THE PMCMC ##########################

# note: log target must do link functions, not the vital rate functions
#       themselves. This for the SDs of parameters which must be +ve

# Functionality for shift to be a par or a constant?

# Setup the prior distributions:
linkdbeta <- function(x, shape1, shape2) dbeta(plogis(x), shape1, shape2)
priorsIPM <- c(rep('dnorm', 11), rep('linkdbeta', 3))

# setup their parameters (note, the names are not required, they just help us
# remember to put things in the right order):
tightPriorsIPM <- list(
  # Survival function intercept:
  sFuncInt=list(mean=-9.65, sd=0.1),
  # Survival function gradient:
  sFuncGra=list(mean=3.77, sd=0.1),
  # Growth function intercept:
  gFuncInt=list(mean=1.41, sd=0.1),
  # Growth function gradient:
  gFuncGra=list(mean=0.56, sd=0.1),
  # Growth function sd:
  gFuncSD=list(mean=log(0.08), sd=1),
  # Reproduction function intercept:
  rFuncInt=list(mean=-7.23, sd=0.1),
  # Reproduction function gradient:
  rFuncInt=list(mean=2.6, sd=0.1),
  # Number of offpsring per litter;
  offNum=list(mean=1, sd=0.1),
  # Childsize intercept:
  osFuncInt=list(mean=0.36, sd=0.1),
  # Childsize gradient:
  osFuncGra=list(mean=0.71, sd=0.1),
  # Childsize sd:
  osFuncSD=list(mean=log(0.16), sd=0.1),
  # Child survival:
  Schild=list(shape1=873, shape2=127),
  # Shift:
  shift=list(shape1=300, shape2=300),
  # The probability of detection of an animal:
  obsProb=list(shape1=800, shape2=200)
)


dispersePriors <- list(
  # Survival function intercept:
  sFuncInt=list(mean=-9.65, sd=1),
  # Survival function gradient:
  sFuncGra=list(mean=3.77, sd=1),
  # Growth function intercept:
  gFuncInt=list(mean=1.41, sd=1),
  # Growth function gradient:
  gFuncGra=list(mean=0.56, sd=1),
  # Growth function sd:
  gFuncSD=list(mean=log(0.08), sd=1),
  # Reproduction function intercept:
  rFuncInt=list(mean=-7.23, sd=1),
  # Reproduction function gradient:
  rFuncInt=list(mean=2.6, sd=1),
  # Number of offpsring per litter;
  offNum=list(mean=1, sd=1),
  # Childsize intercept:
  osFuncInt=list(mean=0.36, sd=1),
  # Childsize gradient:
  osFuncGra=list(mean=0.71, sd=1),
  # Childsize sd:
  osFuncSD=list(mean=log(0.16), sd=1),
  # Child survival:
  Schild=list(shape1=87.3, shape2=12.7),
  # Shift:
  shift=list(shape1=30, shape2=30),
  # The probability of detection of an animal:
  obsProb=list(shape1=80, shape2=20)
)

dispersePriors2 <- list(
  # Survival function intercept:
  sFuncInt=list(mean=-9.65, sd=2),
  # Survival function gradient:
  sFuncGra=list(mean=3.77, sd=2),
  # Growth function intercept:
  gFuncInt=list(mean=1.41, sd=2),
  # Growth function gradient:
  gFuncGra=list(mean=0.56, sd=2),
  # Growth function sd:
  gFuncSD=list(mean=log(0.08), sd=2),
  # Reproduction function intercept:
  rFuncInt=list(mean=-7.23, sd=2),
  # Reproduction function gradient:
  rFuncInt=list(mean=2.6, sd=2),
  # Number of offpsring per litter;
  offNum=list(mean=1, sd=2),
  # Childsize intercept:
  osFuncInt=list(mean=0.36, sd=2),
  # Childsize gradient:
  osFuncGra=list(mean=0.71, sd=2),
  # Childsize sd:
  osFuncSD=list(mean=log(0.16), sd=2),
  # Child survival:
  Schild=list(shape1=87.3, shape2=12.7),
  # Shift:
  shift=list(shape1=30, shape2=30),
  # The probability of detection of an animal:
  obsProb=list(shape1=80, shape2=20)
)

skeleton = list(
  survPars = rep(NA, 2),
  growthPars = rep(NA, 3),
  reprPars = rep(NA, 2),
  offNumPars = NA,
  offSizePars = rep(NA, 3),
  Schild = NA, 
  shift = NA,
  obsProbPar = NA
)

IPMLTP <- list(
  priorFunc = 'evalPriors',
  priors = priorsIPM,
  priorPars = tightPriorsIPM,
  skeleton = skeleton,
  survFunc = 'linLogit',
  growthSamp = 'sampleNorm',
  reprFunc = 'linLogit',
  offNumSamp = 'returnConstant',
  offSizeSamp = 'sampleNorm',
  oneSex = TRUE,
  mu = 'multnomMu',
  muPar = c(trueStart, list(muPs)),
  b = 100, 
  Y = sizeObservations,
  obsProb = 'detectionNumObs',
  breaks = breaks
)

IPMLTPD <- IPMLTP
IPMLTPD$priorPars <- dispersePriors

startValues <- list(
  survPars = c(-9.65, 3.77),
  growthPars = c(1.41, 0.56, log(0.08)),
  survPars = c(-7.23, 2.6),
  offNumPars = 1,
  offSizePars = c(0.36, 0.71, log(0.16)),
  Schild = qlogis(0.873), 
  shift = qlogis(0.44),
  obsProbPar = qlogis(0.8)
)

startValues %<>% unlist
SVjitter <- startValues + rnorm(length(startValues), 0, 0.)

# NOTE: Using this below line of code is a good way to get a general idea of how
#       many particles we require to obtain an estimate of the likelihood which
#       doesn't vary too much:

# Test the specification on the log posterior:
logTargetIPM(startValues, logTargetPars = IPMLTP)
logTargetIPM(SVjitter, logTargetPars = IPMLTP)

# On the disperse priors:
logTargetIPM(startValues, logTargetPars = IPMLTPD)
logTargetIPM(SVjitter, logTargetPars = IPMLTPD)

# Objects which the cluster needs to know about to run the sampling:
clNeeds = c('logTargetIPM', '%<>%', '%>%', 'sampleNorm', 'returnConstant',
            'detectionNumObs', 'evalPriors', 'breaks', 'sampleStateIPM',
            'particleFilter', 'multnomMu', 'linLogit', 'vectorisedSamplerIPM',
            'linkdbeta', 'rmvnorm') 

# Now run a chain:
set.seed(10101)
ungulateChain <- pMH(proposal = multvarNormProp, uFunc = multvarPropUpdate,
                     propPars = diag(length(startValues))/100,
                     lTarg = logTargetIPM, lTargPars = IPMLTPD,
                     x0 = SVjitter, itermax=10000, prntPars = T,
                     clNeeds = clNeeds, RcppFile = "RCPP2.cpp")

# Remove the burn in from the chain, and thin it:
ungulateThinned <- thinMCMC(ungulateChain, alpha = 0.1)

# Plots for a thinned chain:
plot(ungulateThinned, cols = c(2,3), useX11 = T)


############################ A FLAT PRIOR ANALYSIS #############################

priorsFlat <- rep("dunif", 14)

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
  # Shift:
  shift=list(min=-50, max=50),
  # The probability of detection of an animal:
  obsProb=list(min=-50, max=50)
)

plotNames <- c("log(pi)", "Si", "Sg", "Gi", "Gg", "log(Gs)", "Ri", "Rg", "ON",
               "OSi", "OSg", "log(OSs)", "qlogis(Schild)", "qlogis(Shift)",
               "qlogis(ObsP)")

IPMLTP.FLAT <- IPMLTP
IPMLTP.FLAT$priors <- priorsFlat
IPMLTP.FLAT$priorPars <- flatPriors
saveRDS(IPMLTP.FLAT, "IPMLTP.FLAT")

### Note : this loop is iterated a few times to sequentially obtain a 
#          correlation structure for the proposal which converges well enough.
#          The starting sigma was diag(length(startValues))/100, and the number 
#          of times the process was repeated was: 1
set.seed(10101)
ungulateChainFlat <- pMH(proposal = multvarNormProp, uFunc = NULL,
                     propPars = startingSigma,
                     lTarg = logTargetIPM, lTargPars = IPMLTP.FLAT,
                     x0 = SVjitter, itermax=15000, prntPars = T,
                     clNeeds = clNeeds, packages = "dissPackage2")


saveRDS(ungulateChainFlat, "ungChainFlat")
ungulateChainFlat <- readRDS("ungChainFlat")

# Create very high resolution files with all the plots:
plot(ungulateChainFlat, filePath="", cols=2:15, width=8.3*3, height=11.7*3,
     cex=2, names=plotNames)

# Code used to update correlation structure of proposal :
# saveRDS(ungulateChainFlat[[1]][,-1] %>% cov, "flatPriorSigma")
# startingSigma <- readRDS("flatPriorSigma")


thinnedFlat <- thinMCMC(ungulateChainFlat, removeBI = 0.2, alpha = 0.1)

# high resolution plots:
plot(thinnedFlat, filePath="", cols=2:15, width=8.3*3, height=11.7*3,
     cex=2, names=plotNames)

############ RUN AN UNGULATE PMCMCM WITH THE MULTINOMIAL OBS MODEL #############

# load some stored arguments which avoid keeping old code around:
IPMLTP.FLAT <- readRDS("IPMLTP.FLAT")
startValues <- readRDS("startValues")
SVjitter <- readRDS("SVjitter")
clNeeds <- readRDS("clNeeds")

plotNames <- c("log(pi)", "Si", "Sg", "Gi", "Gg", "log(Gs)", "Ri", "Rg", "ON",
               "OSi", "OSg", "log(OSs)", "qlogis(Schild)", "qlogis(Shift)",
               "qlogis(ObsP)")

# Modify the log target parameters such that we now specify the multinomial 
# observation model:
IPMLTPM <- IPMLTP.FLAT
IPMLTPM$obsProb <- 'multinomialObs'

# To avoid copying and pasting large chunks of the same code:
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

# Give us some new data which is more appropriate:
set.seed(102010)
simmedData <- do.call(simulateIBM, simPars)

# Generate the observed count data from the IBM simulations:
max.cens <- simmedData$census.number %>% max
pmcmcBreaks <- seq(1.5, 3.55, l=21)
pmcmcObs <- getSizeDistns(simmedData, pmcmcBreaks)[,(max.cens-5):max.cens]

# Generate a different data set with the same parameters to get a prior 
# on the initial size distribution. We do this in a loop so that we get a 
# better approximation of the true probability of being in a given size set:

priorLoopN <- 100
priorCounts <- rep(0, length(pmcmcBreaks)-1)

set.seed(201020)
for (i in 1:priorLoopN){
  # generate the data:
  priorData <- do.call(simulateIBM, simPars)  
  max.prior <- priorData$census.number %>% max
  # extract the counts in each size class for the first census:
  loop.counts <- subset(priorData, priorData$census.number==(max.prior-5)) %>%
    `$`(size) %>% na.omit %>% vectorToCounts(breaks = pmcmcBreaks)
  # update the previous vector:
  priorCounts %<>% `+`(loop.counts)
}

priorProbs <- countsToProbs(priorCounts)
priorStartSize <- (sum(priorCounts)/priorLoopN) %>% `*`(1.2) %>% ceiling

# Just check some of the values this returns, to see if it's plausible:
multinomialObs(pmcmcObs[,1], multnomMu(priorStartSize, priorProbs), 0.8)

# Setup the rest of our parameter vector:
IPMLTPM$Y <- pmcmcObs
IPMLTPM$breaks <- pmcmcBreaks
IPMLTPM$muPar <- list(priorStartSize, priorProbs)
IPMLTPM$b <-150

# Check some log lik values
logTargetIPM(startValues, logTargetPars = IPMLTPM)
logTargetIPM(SVjitter, logTargetPars = IPMLTPM)

set.seed(10101)
A <- Sys.time()
ungulateChainMult <- pMH(proposal = multvarNormProp, uFunc = NULL,
                         propPars = readRDS("ungMultSigma"),
                         lTarg = logTargetIPM, lTargPars = IPMLTPM,
                         x0 = SVjitter, itermax=5000, prntPars = T,
                         clNeeds = clNeeds, packages = "dissPackage3")
print(Sys.time()-A)

# high resolution plots:
plot(ungulateChainMult, filePath="", cols=2:15, width=8.3*3, height=11.7*3,
     cex=2, names=plotNames)

saveRDS(cov(ungulateChainMult[[1]][,-1]), "ungMultSigma")

####################### TRY WITH THE POISSON MODEL INSTEAD #####################

IPMLTPP <- IPMLTPM
IPMLTPP$obsProb <- 'poissonObs'

# Check some log lik values
logTargetIPM(startValues, logTargetPars = IPMLTPP)
logTargetIPM(SVjitter, logTargetPars = IPMLTPP)


set.seed(10101)
A <- Sys.time()
ungulateChainMult <- pMH(proposal = multvarNormProp, uFunc = NULL,
                         propPars = readRDS("ungMultSigma2"),
                         lTarg = logTargetIPM, lTargPars = IPMLTPM,
                         x0 = SVjitter, itermax=2000, prntPars = T,
                         clNeeds = clNeeds, packages = "dissPackage3")
print(Sys.time()-A)

# high resolution plots:
plot(ungulateChainMult, filePath="", cols=2:15, width=8.3*3, height=11.7*3,
     cex=2, names=plotNames)

saveRDS(cov(ungulateChainMult[[2]][,-1]), "ungMultSigma2")

########################## BACK TO THE COUNTS MODEL ############################

startValues <- readRDS("startValues")
SVjitter <- readRDS("SVjitter")

# shift will be fixed, so remove it:
startValues %<>% as.list
startValues$shift <- NULL
startValues %<>% unlist
SVjitter %<>% as.list
SVjitter$shift <- NULL
SVjitter %<>% unlist

IPMLTP2 <- IPMLTPM
IPMLTP2$breaks <- seq(1.5, 3.55, l=6)[-2]
breaks <- IPMLTP2$breaks
IPMLTP2$obsProb <- "detectionNumObs"
IPMLTP2$b <- 100
IPMLTP2$shift <- qlogis(0.49)

# Get the prior distribution:
priorLoopN <- 100
priorCounts <- rep(0, length(IPMLTP2$breaks)-1)

set.seed(201020)
for (i in 1:priorLoopN){
  # generate the data:
  priorData <- do.call(simulateIBM, simPars)  
  max.prior <- priorData$census.number %>% max
  # extract the counts in each size class for the first census:
  loop.counts <- subset(priorData, priorData$census.number==(max.prior-5)) %>%
    `$`(size) %>% na.omit %>% vectorToCounts(breaks = IPMLTP2$breaks)
  # update the previous vector:
  priorCounts %<>% `+`(loop.counts)
}

priorProbs <- countsToProbs(priorCounts)
priorStartSize <- (sum(priorCounts)/priorLoopN) %>% `*`(1.2) %>% ceiling
IPMLTP2$muPar <- list(priorStartSize, priorProbs)

# Get the observations:
set.seed(102010)
simmedData <- do.call(simulateIBM, simPars)
max.cens <- simmedData$census.number %>% max
IPMLTP2$Y <- getSizeDistns(simmedData, IPMLTP2$breaks)[,(max.cens-5):max.cens]

# Check some log lik values
logTargetIPM(startValues, logTargetPars = IPMLTP2)
logTargetIPM(SVjitter, logTargetPars = IPMLTP2)


set.seed(10101)
A <- Sys.time()
ungulateChainMult <- pMH(proposal = multvarNormProp, uFunc = NULL,
                         propPars = SVjitter %>% length %>% diag %>% `\`(100),
                         lTarg = logTargetIPM, lTargPars = IPMLTP2,
                         x0 = SVjitter, itermax=1000, prntPars = T,
                         clNeeds = clNeeds, packages = "dissPackage3")
print(Sys.time()-A)

# high resolution plots:
plot(ungulateChainMult, filePath="", cols=2:15, width=8.3*3, height=11.7*3,
     cex=2, names=plotNames)

saveRDS(cov(ungulateChainMult[[1]][,-1]), "flatPriorSigma2.0")
saveRDS(ungulateChainMult, 'bigUngChain')
