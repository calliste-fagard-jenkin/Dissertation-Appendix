# Some of the tables and figures in the writeup were produced in this file 
# rather than in the case studies files, to avoid bloat.

# source all the files to have access to the functions written for the
# dissertation:
source('simulateData.R')
source('piecemealFunctions.R')
source('PlottingFunctions.R')
source('ModelSpecIPM.R')
source('SMC.R')
library(ggplot2)
library(doParallel)
# to convert data/frames to LaTeX tables:
library(xtable) 
library(plot.matrix)

############ FOR TABLE 1; AN EXAMPLE OF DATA SIMULATED FROM THE IBM ############

# Produce the data:
set.seed(102030)
SimmedData <- simulateIBM(n=500, t=1000,
                          # set survival details:
                          survFunc = linLogit, survPars = c(-8, 3.77),
                          # set growth details:
                          growthSamp = sampleDTN,
                          growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
                          # set reproduction details:
                          reprFunc = linLogit, reprPars = c(-7, 2.6),
                          # set offspring number and size distribution details:
                          offNum=1, offSizeSamp = sampleDTN,
                          offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
                          # Child survival probability:
                          Schild=0.873,
                          # set other miscelaneous parameters:
                          Start=3, thresh=5000, OneGend = TRUE, popPrint = F)

# Produce the table:
helper1 <- function(x) replace(x, which(is.na(x)), "NA")
SimmedData %>% round(digits=2) %>% lapply(FUN=as.character) %>%
  lapply(FUN=helper1) %>% data.frame %>% xtable


################ FIGURE XX - EXAMPLE TRAJECTORIES OF THE FILTER ################

set.seed(102030)
# We simulate data from a constant size distribution. Projecting this forwards
# allows us to get a bette restimate of the starting distribution to simulate
# from:
sampleN=500
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
                          Start = sampleStart, thresh=5000, OneGend = TRUE,
                          popPrint = T, verbose=F)

# Uncomment the below code to keep only the last 50 censuses:
simmedData <- with(simmedData, # Keep only the last 50 censuses for stability
                   subset(simmedData, census.number>=(max(census.number)-50)))

# Reajust the census.number labels:
simmedData$census.number <- simmedData$census.number -
  min(simmedData$census.number) + 1

######################### MAKE SOME ILLUSTRATIVE PLOTS  ########################

breaks <- seq(1.5, 3.55, l=21)
sizeObservations <- getSizeDistns(simmedData, breaks=breaks)

# Produce a figure:
plot(log(sizeObservations), ylab="log count", xlab='Census number', main='')

# As a placeholder for the prior, use the distribution of observed sizes
# across all censuses:
muPs <- simmedData$size %>% na.omit %>% vectorToCounts(breaks) %>% countsToProbs


# Make a list of the parameters for the particle filter:
ungulatePars <- list(
  survFunc = linLogit, survPars = c(-9.65, 3.77),
  growthSamp = sampleDTN, growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
  reprFunc = linLogit, reprPars = c(-7.23, 2.6), 
  offSizeSamp = sampleDTN, offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
  offNumSamp = returnConstant, offNumPars = 1,
  Schild=qlogis(0.873), oneSex = TRUE, breaks = breaks, shift=qlogis(0.44))

# Simulate a single time step:
do.call(sampleStateIPM, c(list(multnomMu(500, muPs)), ungulatePars, verbose=T))

# Simulate multiple time steps:
obsPopSizes <- getSizeDistns(simmedData, breaks) %>% apply(2, sum)
popSizeAtStart <- subset(simmedData, simmedData$census.number==1)$size %>%
  na.omit %>% length

# Look at the growth rate of the population:
plot(obsPopSizes, type='l', col=4, lwd=3, ylim = c(400, 5000),
     xlab="Census Number", ylab="Population Size")

maxCN <- max(simmedData$census.number)

simSize <- 200
popSizeMatrix <- matrix(NA, nrow=simSize, ncol=length(obsPopSizes)) 

# Work out simulated trajectories:
for (i in 1:simSize){
  popSizeMatrix[i,] <- projectStateSpace(sampleStateIPM, ungulatePars,
      multnomMu(popSizeAtStart, muPs), maxCN-1) %>% apply(2, sum)
}

# Add them to the plot:
meanPopSize <- popSizeMatrix %>% apply(2, mean)
lowerQuant <- popSizeMatrix %>% apply(2, quantile, probs=0.025)
upperQuant <- popSizeMatrix %>% apply(2, quantile, probs=0.975)
lines(meanPopSize, type='l', col=5, lwd=2)
lines(lowerQuant, type='l', col=5)
lines(upperQuant, type='l', col=5)
legend(x='topleft', legend=c("True population size",
      "Simulated Trajectory Mean and 2.5% quantiles"),fill=c(4,5))

########################### COMPARING GRAIN SIZE ###############################
repeats <- 100
simTime <- 10
sizeSet <- c(5, 10, 25, 50, 200)
results <- matrix(NA, nrow = length(sizeSet) + 1, ncol = simTime)

set.seed(101010)

simIBM <- simulateIBM(n=500, t=simTime*10,
                      # set survival details:
                      survFunc = linLogit, survPars = c(-9.65, 3.77),
                      # set growth details:
                      growthSamp = sampleDTN,
                      growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
                      # set reproduction details:
                      reprFunc = linLogit, reprPars = c(-7.23, 2.6),
                      # set offspring number and size distribution deets:
                      offNum=1, offSizeSamp = sampleDTN,
                      offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
                      # Child survival probability:
                      Schild=0.873,
                      # set other miscelaneous parameters:
                      Start=readRDS("startSizes") %>%
                        sample(500, replace=T), thresh=50000,
                      OneGend = TRUE, popPrint = F)

sizeDists <- getSizeDistns(simIBM, breaks)
sizeDists <- sizeDists[,(ncol(sizeDists)-simTime+1):ncol(sizeDists)]
results[1,] <- apply(sizeDists, 2, sum)

for (j in 1:length(sizeSet)){
  resultsj <- matrix(NA, nrow=repeats, ncol=simTime)
  for (i in 1:repeats){
    breaks2 <- seq(1.5, 3.55, l=sizeSet[j]+1)
  
    simPars <- list(
      survFunc = linLogit, survPars = c(-9.65, 3.77),
      growthSamp = sampleDTN, growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
      reprFunc = linLogit, reprPars = c(-7.23, 2.6), 
      offSizeSamp = sampleDTN, offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
      offNumSamp = returnConstant, offNumPars = 1,
      Schild=qlogis(0.873), oneSex = T, breaks = breaks2, shift=qlogis(0.445))
  
    startSize <- results[1, 1]
    sizeDists <- getSizeDistns(simIBM, breaks2)
    sizeDists <- sizeDists[,(ncol(sizeDists)-simTime+1):ncol(sizeDists)]
    startDistn <- sizeDists[,1]/sum(sizeDists[,1])
    resultsj[i,] <- projectStateSpace(sampleStateIPM, simPars,
      multnomMu(startSize, startDistn), simTime - 1) %>% apply(2, sum)
  }
  
  results[j+1, ] <- apply(resultsj, 2, mean)
}

plot(results[2,], type='l', col=adjustcolor(2, alpha=0.5), lwd=2,
     ylab='Population size', lty=2,
     xlab="Census number", ylim = c(4500, 6600))

for (i in 2:length(sizeSet)){
  lines(results[i+1,], col=adjustcolor(col=i+1, alpha=0.5), lwd=2, lty=2)
}

legend('topleft', fill = adjustcolor(col = 2:(length(sizeSet)+1), alpha = 0.5),
       legend = as.character(sizeSet), bty = 'n', horiz = T,
       title = "Number of Size Classes")


####### MAKE A PLOT WITH SIMULATIONS FROM NORMAL AND STATE SPACE SAMPLER #######

simSize <- 50
simTime <- 10
IBMresults <- matrix(NA, nrow=simTime, ncol=simSize)
SSMresults <- matrix(NA, nrow=simTime, ncol=simSize)
SSMresults2 <- matrix(NA, nrow=simTime, ncol=simSize)

set.seed(101010)
for (i in 1:simSize){
  # produce the IBM simulation:
  simIBM <- simulateIBM(n=500, t=simTime*10,
                            # set survival details:
                            survFunc = linLogit, survPars = c(-9.65, 3.77),
                            # set growth details:
                            growthSamp = sampleDTN,
                            growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
                            # set reproduction details:
                            reprFunc = linLogit, reprPars = c(-7.23, 2.6),
                            # set offspring number and size distribution deets:
                            offNum=1, offSizeSamp = sampleDTN,
                            offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
                            # Child survival probability:
                            Schild=0.873,
                            # set other miscelaneous parameters:
                            Start=readRDS("startSizes") %>%
                            sample(500, replace=T), thresh=50000,
                            OneGend = TRUE, popPrint = F)
  
  # extract the starting size distribution:
  sizeDists <- getSizeDistns(simIBM, breaks)
  sizeDists <- sizeDists[,(ncol(sizeDists)-simTime+1):ncol(sizeDists)]
  IBMresults[,i] <- apply(sizeDists, 2, sum)
  startSize <- IBMresults[1, i]
  startDistn <- sizeDists[,1]/sum(sizeDists[,1])
  
  # Set the SSM parameters:
  simPars <- list(
    survFunc = linLogit, survPars = c(-9.65, 3.77),
    growthSamp = sampleDTN, growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
    reprFunc = linLogit, reprPars = c(-7.23, 2.6), 
    offSizeSamp = sampleDTN, offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
    offNumSamp = returnConstant, offNumPars = 1,
    Schild=qlogis(0.873), oneSex = TRUE, breaks = breaks, shift=qlogis(0.445))
  
  
  # simulate the results with the right shift:
  SSMresults[,i] <- projectStateSpace(sampleStateIPM, simPars,
        multnomMu(startSize, startDistn), simTime - 1) %>% apply(2, sum)
  
  # simulate the results with an incorrect shift of 0.5 (qlogis(0)=0.5):
  simPars$shift <- 0
  SSMresults2[,i] <- projectStateSpace(sampleStateIPM, simPars,
     multnomMu(startSize, startDistn), simTime - 1) %>% apply(2, sum)
}

# Calculate the 95% lines for the IBM data:
meanIBM <- apply(IBMresults, 1, mean)
lowerIBM <- apply(IBMresults, 1, quantile, probs=0.025)
upperIBM <- apply(IBMresults, 1, quantile, probs=0.975)

# Calculate the 95% lines for the SSM data with the correct shift:
meanSSM <- apply(SSMresults, 1, mean)
lowerSSM <- apply(SSMresults, 1, quantile, probs=0.025)
upperSSM <- apply(SSMresults, 1, quantile, probs=0.975)

# Calculate the 95% lines for the SSM data with the incorrect shift:
meanSSM2 <- apply(SSMresults2, 1, mean)
lowerSSM2 <- apply(SSMresults2, 1, quantile, probs=0.025)
upperSSM2 <- apply(SSMresults2, 1, quantile, probs=0.975)

# Set some colours to reuse:
col4 <- adjustcolor(4, alpha=0.8)
col5 <- adjustcolor(5, alpha=0.8)


# First figure (with correct shift):
plot(meanIBM, type='l', lwd=2, col=col4, ylim=c(2500, 7000),
     ylab='Population size', xlab='Census number')

lines(lowerIBM, col=col4, lty=2) ; lines(upperIBM, col=col4, lty=2)
lines(lowerSSM, col=col5, lty=2) ; lines(upperSSM, col=col5, lty=2)
lines(meanSSM, col=col5, lwd=2)
legend('topleft', legend=c("IBM simulations", "SSM simulations"),
       fill=c(col4,col5), bty='n')

# Second figure (with incorrect shift):
plot(meanIBM, type='l', lwd=2, col=col4, ylim=c(2500, 7100),
     ylab='Population size', xlab='Census number')

lines(lowerIBM, col=col4, lty=2) ; lines(upperIBM, col=col4, lty=2)
lines(lowerSSM2, col=col5, lty=2) ; lines(upperSSM2, col=col5, lty=2)
lines(meanSSM2, col=col5, lwd=2)
legend('topleft', legend=c("IBM simulations", "SSM simulations"),
       fill=c(col4,col5), bty='n')

############## STACKED AREA GRAPH FOR OVIS ARIES PMCMC DATA ####################
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

# generate the data again, so that we can plot it for the whole time series:
set.seed(102010)
simmedData <- do.call(simulateIBM, simPars)

# make the plot:
getSizeDistns(simmedData, seq(1.5, 3.55, l=6)) %>% 
  plotSizeCounts(breaks = seq(1.5, 3.55, l=6))

# make the table with 5 classes:
colnames(ungData) <- paste("census", 1:6)
rownames(ungData) <- paste("size class", 1:5)
xtable(ungData, digits = 0)

# make the table with 4 classes:
ungData[1:2,] %>% apply(2, sum) %>% rbind(ungData[3:5,]) %>% xtable(digits = 0)


############### FINDING OUT SHIFT FOR UNGULATE PMCMC DATA ######################

breaks <- seq(1.5, 3.55, l=6)[-2]

simSize <- 20
simTime <- 10
IBMresults <- matrix(NA, nrow=simTime, ncol=simSize)
SSMresults <- matrix(NA, nrow=simTime, ncol=simSize)

simParsPlot <-  simPars
simParsPlot$n <- 200
simParsPlot$t <- 140
simParsPlot$thresh <- NULL
simParsPlot$Start <-  3

set.seed(101010)
for (i in 1:simSize){
  # produce the IBM simulation:
  simIBM <- do.call(simulateIBM, simParsPlot)
  
  # extract the starting size distribution:
  sizeDists <- getSizeDistns(simIBM, breaks)
  sizeDists <- sizeDists[,(ncol(sizeDists)-simTime+1):ncol(sizeDists)]
  IBMresults[,i] <- apply(sizeDists, 2, sum)
  startSize <- IBMresults[1, i]
  startDistn <- sizeDists[,1]/sum(sizeDists[,1])
  
  # Set the SSM parameters:
  simPars2 <- list(
    survFunc = linLogit, survPars = c(-9.65, 3.77),
    growthSamp = sampleDTN, growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
    reprFunc = linLogit, reprPars = c(-7.23, 2.6), 
    offSizeSamp = sampleDTN, offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
    offNumSamp = returnConstant, offNumPars = 1,
    Schild=qlogis(0.873), oneSex = TRUE, breaks = breaks, shift=qlogis(0.50))
  
  # simulate the results with the right shift:
  SSMresults[,i] <- projectStateSpace(sampleStateIPM, simPars2,
     multnomMu(startSize, startDistn), simTime - 1) %>% apply(2, sum)
}

# Calculate the 95% lines for the IBM data:
meanIBM <- apply(IBMresults, 1, mean)
lowerIBM <- apply(IBMresults, 1, quantile, probs=0.025)
upperIBM <- apply(IBMresults, 1, quantile, probs=0.975)

# Calculate the 95% lines for the SSM data with the correct shift:
meanSSM <- apply(SSMresults, 1, mean)
lowerSSM <- apply(SSMresults, 1, quantile, probs=0.025)
upperSSM <- apply(SSMresults, 1, quantile, probs=0.975)

# Set some colours to reuse:
col4 <- adjustcolor(4, alpha=0.8)
col5 <- adjustcolor(5, alpha=0.8)

# First figure (with correct shift):
plot(meanIBM, type='l', lwd=2, col=col4, ylim=c(500, 11000),
     ylab='Population size', xlab='Census number')

lines(lowerIBM, col=col4, lty=2) ; lines(upperIBM, col=col4, lty=2)
lines(lowerSSM, col=col5, lty=2) ; lines(upperSSM, col=col5, lty=2)
lines(meanSSM, col=col5, lwd=2)
legend('topleft', legend=c("IBM simulations", "SSM simulations"),
       fill=c(col4,col5), bty='n')

############### MAKE PLOT FOR LOG POSTERIOR SD FOR DIFFERENT B #################

IPMLTP <- readRDS("IPMLTP")

startValues <- list(
  survPars = c(-9.65, 3.77),
  growthPars = c(1.41, 0.56, log(0.08)),
  survPars = c(-7.23, 2.6),
  offNumPars = 1,
  offSizePars = c(0.36, 0.71, log(0.16)),
  Schild = qlogis(0.873), 
  obsProbPar = 15
)

startValues %<>% unlist

# setup for the loop:
bVals <- c(10, 25, 50, 100, 150, 200, 250, 500)
simNum <- 50
results <- matrix(NA, nrow = length(bVals), ncol = simNum)
IPMLTP.graph <- IPMLTP

# running the simulation for estimate MC variance:
for (i in 1:length(bVals)){
  IPMLTP.graph$b <- bVals[i]
  for (j in 1:simNum){
    results[i, j] <- logTargetIPM(startValues, logTargetPars = IPMLTP.graph,
                                  returnLL = T)
  }
  cat("Loop for b =",bVals[i],"has been completed\n")
}

results %>% apply(1, sd) %>% plot(x = bVals, xlab = "Number of particles", 
                                  ylab = "Standard deviation of log likelihood",
                                  col = "purple", pch = 16)

lines(x = bVals, y = results %>% apply(1, sd), col = "purple", lwd = 2)

#################### ESTIMATE SD OF LOG LIKELIHOOD IN PARALLEL #################
IPMLTP.copy <- IPMLTP
IPMLTP.copy$b <- 500
llSD <-  rep(NA, 100)
cl<-makeCluster(4)
registerDoParallel(cl)
clusterExport(cl, c("%<>%", "%>%", "detectionNumObs", "evalPriors", "multnomMu",
                    "linLogit", "sampleNorm", "returnConstant"))

# Do the parallel loop:
ls <- foreach(i = 1:length(llSD), .packages = 'dissPackage3') %dopar%{
  # then run the chain on the cluster:
  to.ls <- logTargetIPM(startValues, logTargetPars = IPMLTP.copy, returnLL=T)
}

try(stopCluster(cl), silent = T)
ls %>% unlist %>% sd

############## GET EFFECTIVE SAMPLE SIZE OF PARTICLE FILTER ####################
simNum <-  500
IPMLTP.copy <- IPMLTP
IPMLTP.copy$b <- 100
cl<-makeCluster(4)
registerDoParallel(cl)
clusterEvalQ(cl, library("dissPackage3"))
clusterExport(cl, c("%<>%", "%>%", "detectionNumObs", "evalPriors", "multnomMu",
                    "linLogit", "sampleNorm", "returnConstant"))

ls <-  foreach(i = 1:simNum) %dopar%{
  ll <- logTargetIPM(startValues, logTargetPars = IPMLTP.copy, returnW = T) %>%
    apply(2, countsToProbs)
  to.ls <- ll %>% `^`(2) %>% apply(2, sum) %>% `^`(-1) %>% `[`(-1) %>% mean
}

ls %>% unlist %>% mean
try(stopCluster(cl), silent = T)

##################### SHEEP DATA DOUCET SD OF LOG LIK ##########################

sheepIPMLTP <- readRDS('sheepIPMLTP')
sheepMLEs <- readRDS('sheepMLEs')

IPMLTP.copy <- sheepIPMLTP
IPMLTP.copy$b <- 100
llSD <-  rep(NA, 100)
cl<-makeCluster(4)
registerDoParallel(cl)
clusterExport(cl, c("%<>%", "%>%", "detectionNumObs", "evalPriors", "multnomMu",
                    "linLogit", "sampleNorm", "returnConstant"))

# Do the parallel loop:
ls <- foreach(i = 1:length(llSD), .packages = 'dissPackage3') %dopar%{
  # then run the chain on the cluster:
  to.ls <- logTargetIPM(sheepMLEs$par, logTargetPars = IPMLTP.copy, returnLL=T)
}

try(stopCluster(cl), silent = T)
ls %>% unlist %>% sd