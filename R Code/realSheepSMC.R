source('AdaptivePMCMC.R')
source('simulateData.R')
source('ModelSpecIPM.R')
source('piecemealFunctions.R')
source('SMC.R')
library(dissPackage3)
library(MASS)
library(mcmcse)

SHEEP <- read.csv("SHEEP.csv", header = 2)

# How many years of data:
SHEEP$sheep.yr %>% unique %>% length

# Remove all dead sheep:
SHEEP.ALIVE <- subset(SHEEP, !is.na(SHEEP$wtt1))
SHEEP.ALIVE.O <- SHEEP.ALIVE[order(SHEEP.ALIVE$sheep.yr),]

# Add some handy columns to use a function that exists in another file:
SHEEP.ALIVE$census.number <- SHEEP.ALIVE$sheep.yr-min(SHEEP.ALIVE$sheep.yr)+1
SHEEP.ALIVE$size <- SHEEP.ALIVE$wtt1
SHEEP$census.number <- SHEEP$sheep.yr-min(SHEEP$sheep.yr)+1
SHEEP$size <- SHEEP$wtt1

# Get breakpoints, count data, and growth rate per census:
brks <- c(0, 16, 20, 24, 28, 40)
SHEEP.COUNTS <- getSizeDistns(SHEEP, brks)
SHEEP.COUNTS %>% apply(2, sum)
SHEEP.COUNTS %>% apply(2, sum) %>% divv

popSizes <-  rep(NA, ncol(SHEEP.COUNTS))
for (i in 1:ncol(SHEEP.COUNTS)){
  popSizes[i] <- subset(SHEEP, SHEEP$sheep.yr==(1985+i))[1,]$pop.size
}

obsPopSizes <- rep(NA, ncol(SHEEP.COUNTS))
for (i in 1:ncol(SHEEP.COUNTS)){
  obsPopSizes[i] <- subset(SHEEP, SHEEP$sheep.yr==(1985+i) &
                             !is.na(SHEEP$surv)) %>% nrow
}

popCount <- rep(NA, ncol(SHEEP.COUNTS))
for (i in 1:ncol(SHEEP.COUNTS)){
  popCount[i] <- subset(SHEEP, SHEEP$census.number==i &
                             !is.na(SHEEP$size)) %>% nrow
}


# Try and find out what proportion of the population of females gets detected:
# We make the assumption that detectability of females is 1 (and have confirmed
# with experts that this is reasonable). So in fact, we know the true population
# size at each time step, and the state space model is less appropriate.
# However, we would like to see what the methodology produces in this context.
detectedLiveFemales <- rep(NA, ncol(SHEEP.COUNTS))
for (i in 1:ncol(SHEEP.COUNTS)){
  detectedLiveFemales[i] <- subset(SHEEP, SHEEP$sheep.yr==(1985+i) &
                             SHEEP$surv==1) %>% nrow
}


# sheep counts:
plot(SHEEP.COUNTS %>% apply(2, sum), pch = 16, col = "purple", x = 1986:1996,
     ylab = "population size", xlab = "Sheep year", ylim = c(20, 280))
lines(SHEEP.COUNTS %>% apply(2, sum), col = "purple", x = 1986:1996)
# pop.size column:
lines(detectedLiveFemales, x = 1986:1996, pch = 16, type='p', col = 'Blue')
lines(detectedLiveFemales, x = 1986:1996, col = 'Blue')
# add key:
legend('topleft',legend=c("Captured and weighed surviving females",
                          "Surviving females"), bty = 'n',
       fill = c("purple", "blue"))

# Get the observed capture prob at each time step, fit a beta dist:
captureProb <- apply(SHEEP.COUNTS, 2, sum)/detectedLiveFemales
pars <- fitdistr(captureProb,dbeta,start=list(shape1=1,shape2=1))

# plot the fitted dist over the data observed capture probs:
x <-  seq(0, 1, l=3000)
y <- do.call(dbeta, c(list(as.numeric(x)), as.list(pars$estimate)))
hist(captureProb, freq = F, xlab = "Capture probability", main = "")
lines(x, y, col = "blue", lwd = 2)
points(captureProb, y=rep(0, length(captureProb)), col = "blue", pch=16)

# Get the pars of the initial size dist:
sheepStartSize <- detectedLiveFemales[1]
sheepMuPar <- SHEEP.COUNTS[,1] %>% countsToProbs

########################### SETUP THE PMCMC ####################################

# A function that evaluates the prior for functions that have a logis link and
# a beta prior:
linkdbeta <- function(x, shape1, shape2) dbeta(plogis(x), shape1, shape2)

# Setup the prior names:
priorsIPM <- c(rep("dunif", 12), rep("linkdbeta", 2))

# Setup the prior hyper par values:
sheepPriors <- list(
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
  obsProb=list(shape1=pars$estimate[1], shape2=pars$estimate[2]),
  # shift:
  shift=list(shape1=5, shape2=5)
)

skeleton = list(
  survPars = rep(NA, 2),
  growthPars = rep(NA, 3),
  reprPars = rep(NA, 2),
  offNumPars = NA,
  offSizePars = rep(NA, 3),
  Schild = NA,
  obsProbPar = NA,
  shift = NA
)

IPMLTP <- list(
  priorFunc = 'evalPriors',
  priors = priorsIPM,
  priorPars = sheepPriors,
  skeleton = skeleton,
  survFunc = 'linLogit',
  growthSamp = 'sampleNorm',
  reprFunc = 'linLogit',
  offNumSamp = 'returnConstant',
  offSizeSamp = 'sampleNorm',
  oneSex = TRUE,
  mu = 'multnomMu',
  muPar = c(sheepStartSize, list(sheepMuPar)),
  b = 100, 
  Y = SHEEP.COUNTS,
  obsProb = 'detectionNumObs',
  breaks = brks
)

startValues <- list(
  survPars = c(5.26, 2.84),
  growthPars = c(8.53, -1.74, -2.3),
  reprPars = c(-3.0, -1.5),
  offNumPars = 1.13,
  offSizePars = c(0.3, -2.2, -0.1),
  Schild = 0.9, 
  obsProbPar = -0.25, # not too close to 50 since this will hurt the chain
  shift = 0.44
)

startValues %<>% unlist
logTargetIPM(startValues, logTargetPars = IPMLTP, returnNeg = T, printProp = F)

IPMLTP2 <- IPMLTP
IPMLTP2$b <- 250
# get the MLEs by using more iterations:
MLEs <- optim(startValues, logTargetIPM, logTargetPars = IPMLTP2, returnNeg = T,
              method = "SANN", control = list(maxit = 1000), printProp = T)

# saveRDS(MLEs, "SheepMLEs2")

# Check some log lik values
logTargetIPM(startValues, logTargetPars = IPMLTP, returnNeg = T, printProp = F)
logTargetIPM(MLEs$par, logTargetPars = IPMLTP, returnNeg = T, printProp = F)

# get the hessian:
A <- Sys.time()
SANN <-  optim(MLEs$par, logTargetIPM, logTargetPars = IPMLTP, returnNeg = T,
               method = "SANN", control = list(maxit = 10), printProp = T,
               hessian = T)
print(Sys.time()-A)

# get the covariance of the proposal:
propCov <- solve(SANN$hessian)
# saveRDS(propCov, "sheepPropCov")


# high resolution plots:
plotNames <- c("log(pi)", "Si", "Sg", "Gi", "Gg", "log(Gs)", "Ri", "Rg", "ON",
               "OSi", "OSg", "log(OSs)", "qlogis(Schild)", "qlogis(ObsP)",
               "qlogis(shift)")

# Objects which the cluster needs to know about to run the sampling:
clNeeds = c('logTargetIPM', '%<>%', '%>%', 'sampleNorm', 'returnConstant',
            'detectionNumObs', 'evalPriors', 'sampleStateIPM', 'particleFilter',
            'multnomMu', 'linLogit', 'vectorisedSamplerIPM', 'rmvnorm',
            'linkdbeta')

# the functions for the backtransform of each parameter:
returnSelf <- function(x) x
links <- c(rep('returnSelf', 4), 'exp', rep('returnSelf', 5), 'exp',
           rep('plogis', 3))

############################## RUN THE PMCMC ###################################
IPMLTP$b <- 100
sheepStart <- readRDS("sheepChain4")[[1]][10000,-1]
sheepChain <- pMH(proposal = multvarNormProp, uFunc = multvarPropUpdate,
                         propPars = diag(length(startValues))/60,
                         #propPars = readRDS("sheepSD8"),
                         #propPars = propCov,
                         cores = 3, nChains = 3,
                         lTarg = logTargetIPM, lTargPars = IPMLTP,
                         x0 = sheepStart, itermax=24000, prntPars = T,
                         clNeeds = clNeeds, packages = "dissPackage3")

getAcceptance(sheepChain, start = 12000)

# produce trace, acf, and conditional marginal plots:
plot(sheepChain, cols=2:15, cex=2, names=plotNames,
     width=8.3*3, height=11.7*3, filePath="sheepPlots")

# thin chain and produce a version without the burn in:
thinnedSheepChain <- thinMCMC(sheepChain, alpha = 0.05, removeBI = 0.6)
sheepChainNoBI <- thinMCMC(sheepChain, alpha = 1, removeBI = 0.6)

# effective sample sizes for the chains:
multiESS(thinnedSheepChain[[2]][,-1]) ; multiESS(thinnedSheepChain[[3]][,-1]) 
multiESS(sheepChainNoBI[[2]][,-1]) ; multiESS(sheepChainNoBI[[3]][,-1]) 

# trace, acf, marginal plots for thinned chain:
plot(thinnedSheepChain, cols=2:15, cex=2, names=plotNames,
     width=8.3*3, height=11.7*3, filePath="sheepPlotsThinned")

# posterior growth rate from the samples:
cluster <- makeCluster(4)
CI <- posteriorGrowthRate(list(thinnedSheepChain[[1]], thinnedSheepChain[[2]]),
                          IPMLTP, growthFunc = "normal", printRate = 10,
                          offSizeFunc = "normal", L = 0, U = 40, m = 500,
                          cluster = cluster, clNeeds = clNeeds)
par(mfrow = c(4,4))
grDiagnostic(sheepChainNoBI)

################# PLOTS ON FIRST 2 PARAMETERS FOR WRITE-UP #####################

# make plots for the first three parameters as figures for the writeup:
chainFirst2 <- rep(NA, 2) %>% as.list
for (j in 2:3) chainFirst2[[j-1]] <- sheepChainNoBI[[j]][,-(4:15)] 
class(chainFirst2) <- 'pMCMCoutput'
F2plotNames <- plotNames[1:3]
plot(chainFirst2, cols = 2:3, width=8.3*3*(0.4), height=11.7*1.5*(0.3),
     cex=1, names=F2plotNames, filePath="chainFirst2_")
grDiagnostic(chainFirst2)

for (j in 2:3) chainFirst2[[j-1]] <- thinnedSheepChain[[j]][,-(4:15)]
class(chainFirst2) <- 'pMCMCoutput'
plot(chainFirst2, cols = 2:3, width=8.3*3*(0.4), height=11.7*1.5*(0.3),
     cex=1, names=F2plotNames, filePath="chainFirst2_")

