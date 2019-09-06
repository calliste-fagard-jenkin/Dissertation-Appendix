source('simulateData.R')
source('piecemealFunctions.R')
source('PlottingFunctions.R')
library(xtable)

# Taken from 2.6 of Ellner, Childs & Rees, this file tries to use various
# pieces of other code to simulate a case study of Ovis Aries, a sheep found
# in Scotland. It was preferred over the case study in chapter 2.5 as the life
# cycle and reproductive behaviour is closer to that of the Guppy population we
# seek to model later.

############################# PRODUCE DATA #####################################
set.seed(102030)
simmedData <- simulateIBM(n=500, t=1000,
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

# Take a sample of the data like the IPM book:
final.census.number <- max(simmedData$census.number)
solveDF <- subset(simmedData, simmedData$census.number==final.census.number)
solveDF <- solveDF[solveDF %>% nrow %>% sample(size=3000) %>% sort,]

# print out mean reproductive rate:
simmedData$reproduced %>% na.omit %>% mean
############### PRODUCE SOME DATA FRAMES FOR EASY MODEL FITTING ################
solveDFg <- with(solveDF, subset(solveDF, !is.na(size) & !is.na(prev.size)))

# remove the observations we can't use for survival fitting:
solveDF2 <- with(solveDF, subset(solveDF, !is.na(prev.size)))

# Make a DF with only the individuals when they were born:
solveDF3 <- subset(solveDF, !is.na(solveDF$parent.size))

# Make a DF with only individuals that had the chance to reproduce:
solveDF4 <-  subset(solveDF, !is.na(solveDF$reproduced))

# Make DF to estimate child survival probability:
solveDF5 <- solveDF[(solveDF$reproduced & !is.na(solveDF$reproduced)) %>% which,]

############## PIECEMEAL ESTIMATION OF DEMOGRAPHIC PARAMETERS ##################

# Estimation without truncation:
growthSol <- lm(size ~ prev.size, data = solveDFg)
survivalSol <- glm(survived ~ prev.size, family = binomial, data = solveDF2)
offspringSizeSol <- lm(size ~ parent.size, data = solveDF3)
reproductionSol <- glm(reproduced ~ size, family = binomial, data=solveDF4)
offspringSurvSol <- glm(off.survived~1, family=binomial, data=solveDF5)

# Estimation with truncation:
growthSolStart <- c(coef(growthSol), summary(growthSol)$sigma)
growthSolTrunc <- optim(growthSolStart, growthNLL, DF=solveDFg, L=1.5, U=3.55,
                        formula=size ~ prev.size)
offspringSizeStart <- c(coef(offspringSizeSol), summary(offspringSizeSol)$sigma)
offspringSizeSolTrunc <- optim(offspringSizeStart, growthNLL, DF=solveDF3,
                           formula=size~parent.size, L=1.5, U=3.55)

# Uncomment below lines to see parameter estimates:
# survivalSol %>% coef
# growthSol %>% coef
# offspringSizeSol %>% coef
# offspringSurvSol %>% coef %>% plogis
# reproductionSol %>% coef

# parameters that don't involve truncated dists:
survival.i <- c(-9.65, coef(survivalSol)[1])
survival.g <- c(3.77, coef(survivalSol)[2])
repro.i <- c(-7.23, coef(reproductionSol)[1])
repro.g <- c(2.6, coef(reproductionSol)[2])
offsurv.p <- c(0.873, offspringSurvSol %>% coef %>% plogis)

# parameters that do involve truncated dists:
growth.i <- c(1.41, coef(growthSol)[1], growthSolTrunc$par[1])
growth.g <- c(1.56, coef(growthSol)[2], growthSolTrunc$par[2])
growth.s <- c(0.08, summary(growthSol)$sigma, exp(growthSolTrunc$par[3]))
offsize.i <- c(0.36, coef(offspringSizeSol)[1], offspringSizeSolTrunc$par[1])
offsize.g <- c(0.71, coef(offspringSizeSol)[2], offspringSizeSolTrunc$par[2])
offsize.s <- c(0.36, summary(offspringSizeSol)$sigma,
               exp(offspringSizeSolTrunc$par[3]))

noTruncDF <- data.frame(survival.i = survival.i, survival.g = survival.g,
                        repro.i = repro.i, repro.g = repro.g,
                        offsurv.p = offsurv.p)

truncDF <- data.frame(growth.i = growth.i, growth.g = growth.g,
                      growth.s = growth.s, offsize.i = offsize.i, 
                      offsize.g = offsize.g, offsize.s = offsize.s)

xtable(noTruncDF, digits = 4)
xtable(truncDF, digits = 4)

########### PRODUCE A PLOT TO VISUALISE THE REPRODUCTION FUNCTION ##############
xsize <- seq(1.5, 3.55, length=1000)
yreprEst <- predict(reproductionSol, list(size=xsize), type='response')
yreprTru <- linLogit(xsize, c(-7.23, 2.6))

# plot the observations:
plot(solveDF4$size, solveDF4$reproduced, col=adjustcolor('black', alpha=0.05),
     pch=16, xlab='log(size) of individual', ylab='Probability of reproduction')

# add the fitted and true curve:
lines(xsize, yreprEst)
lines(xsize, yreprTru, col='blue')

# calculate the percentage error for each point:
percErr <- abs(yreprEst-yreprTru)/yreprTru
lines(xsize, percErr, col='red')
legend('topleft', legend=c('Fitted Curve', 'True Curve', 'Relative Error'),
       fill=c('black','blue', 'red'), cex=0.8, bty='n')

# same plot for survival:
plot(solveDF2$prev.size, solveDF2$survived, col=adjustcolor('black', alpha=0.05),
     pch=16, xlab='log(size) of individual', ylab='Probability of survival')

# true and est curves:
ysurvEst <- predict(survivalSol, list(prev.size=xsize), type='response')
ysurvTru <- linLogit(xsize, c(-9.65, 3.77))
lines(xsize, ysurvEst) ; lines(xsize, ysurvTru, col='blue')

# perc err curve:
percErr <- abs(ysurvEst-ysurvTru)/ysurvTru
lines(xsize, percErr, col='red')

legend('topleft', legend=c('Fitted Curve', 'True Curve', 'Relative Error'),
       fill=c('black','blue', 'red'), cex=0.8, bty='n')

################ MIDPOINT RULE ITERATION MATRIX APPROXIMATION ##################

# Calculate the projection matrix using estimated parameters:
A <- kernelOneVar(m=500,
                  # growth:
                  growthFunc = doublyTruncatedNormal,
                  growthPars = c(coef(growthSol),
                                 growthSol %>% summary %>% `$`(sigma) %>% log,
                                 1.5, 3.55),
                  # survival:
                  survFunc = linLogit,
                  survPars = coef(survivalSol),
                  # reproduction:
                  repFunc = linLogit,
                  repPars = coef(reproductionSol),
                  # offspring number per birth:
                  offNum = 1,
                  # offpspring sizes:
                  offSizeFunc = doublyTruncatedNormal,
                  offSizePars = c(coef(offspringSizeSol), offspringSizeSol %>%
                                 summary %>% `$`(sigma) %>% log, 1.5, 3.55),
                  # child survival rate:
                  childSurv = offspringSurvSol %>% coef %>% plogis,
                  # Other parameters:
                  halfPop = TRUE, L=1.5, U=3.55)

# Get the projection matrix for the true simulated values:
Atrue <- kernelOneVar(m=500,
                      # Growth Function specification:
                      growthFunc = doublyTruncatedNormal,
                      growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
                      # Survival function specification:
                      survFunc = linLogit, survPars = c(-9.65, 3.77),
                      # Reproductive probability function specification:
                      repFunc = linLogit, repPars = c(-7.23, 2.60),
                      # Offspring distribution functions specification:
                      offNum = 1, offSizeFunc = normal,
                      offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
                      childSurv = 0.873,
                      # Specification of other general parameters:
                      halfPop = TRUE, L=1.5, U=3.55)

# Get the projection matrix for the exact functions and values used in the book:
Abook <- kernelOneVar(m=500,
                      # Growth Function specification:
                      growthFunc = normal,
                      growthPars = c(1.41, 0.557, log(0.0799), 1.6, 3.7),
                      # Survival function specification:
                      survFunc = linLogit, survPars = c(-9.65, 3.77),
                      # Reproductive probability function specification:
                      repFunc = linLogit, repPars = c(-7.23, 2.60),
                      # Offspring distribution functions specification:
                      offNum = 1, offSizeFunc = normal,
                      offSizePars = c(0.362, 0.709, log(0.159), 1.6, 3.7),
                      childSurv = plogis(1.93),
                      # Specification of other general parameters:
                      halfPop = TRUE, L=1.6, U=3.7)

A %>% plot(L=1.5, U=3.55, cex=1.1)
Atrue %>% plot(L=1.5, U=3.55, cex=1.1)
# Abook %>% plot(L=1.5, U=3.55)

# Calculate the growth rates:
A %>% eigen %>% `$`(values) %>% `[`(1) %>% Re
Atrue %>% eigen %>% `$`(values) %>% `[`(1) %>% Re
Abook %>% eigen %>% `$`(values) %>% `[`(1) %>% Re

######## INTEGRATION GRAIN SIZE EFFECT ON GROWTH RATE ESTIMATION ###############

mSizes <- seq(20, 1000, length=10)
growthRates <- length(mSizes) %>% rep(x=NA)

for (i in 1:length(mSizes)){
  B <- kernelOneVar(m=mSizes[i],
                    # growth:
                    growthFunc = doublyTruncatedNormal,
                    growthPars = c(coef(growthSol),
                                   growthSol %>% summary %>% `$`(sigma) %>% log,
                                   1.5, 3.55),
                    # survival:
                    survFunc = linLogit,
                    survPars = coef(survivalSol),
                    # reproduction:
                    repFunc = linLogit,
                    repPars = coef(reproductionSol),
                    # offspring number per birth:
                    offNum = 1,
                    # offpspring sizes:
                    offSizeFunc = doublyTruncatedNormal,
                    offSizePars = c(coef(offspringSizeSol), offspringSizeSol %>%
                                      summary %>% `$`(sigma) %>% log, 1.5, 3.55),
                    # child survival rate:
                    childSurv = offspringSurvSol %>% coef %>% plogis,
                    # Other parameters:
                    halfPop = TRUE, L=1.5, U=3.55)
  
  growthRates[i] <- B %>% eigen %>% `$`(values) %>% `[`(1) %>% Re
}

growthRates %>% plot(x=mSizes, type='p', ylab='Approximated growth rate',
                     xlab='Midpoint rule meshpoints', col='purple', pch=16,
                     cex=1, cex.lab=1.2)

############ GET IDEA OF SIZE DISTRIBUTION OF UNGULATES OVER TIME ##############

censusNums <- simmedData$census.number
simmedData %>% subset(censusNums==max(censusNums)) %>% `$`(size) %>% na.omit %>%
  hist