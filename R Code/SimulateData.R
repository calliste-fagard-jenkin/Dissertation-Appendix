library(magrittr) # for the ever useful pipe operator

linLogit <- function(z, par){
  # purpose : Calculates the probability an individual of size x survives to 
  #           the next time step. Uses a logistic form.
  # inputs  : z         - The size of the individual (continuous)
  #           par       - A vector with entries 'intercept' and 'gradient'
  # output  : The density of survival
  
  # extract parameters:
  intercept <- par[1]
  gradient <- par[2]
  # calculate probability of survival:
  s <- intercept + gradient*z
  plogis(s) %>% return
}


sampleRec <- function(n, rate, L, U){
  # purpose : sample sizes from the recruitment function using a rejection 
  #           algorithm.
  # inputs  : n    - The desired number of samples
  #           rate - The rate parameter of the function
  #           L    - The lower bound of the size variable 
  #           U    - The upper bound of the size variable
  #
  # NOTE : log link on rate?
  samples <- rep(NA, n)
  for (i in 1:n){
    sampled.size <- rexp(1, rate)
    while(sampled.size<L | sampled.size>U) sampled.size <- rexp(1, rate)
    samples[i] <- sampled.size
  }
  
  return(samples)
}

sampleDTN <- function(x, pars){
  # purpose : produces samples from a doubly truncated normal distribution 
  #           using a rejection sampling algorithm
  # extract parameter values:
  if (length(pars)!=5) stop('Incorrect number of parameters')
  int <- pars[1] ; gra <- pars[2] ; sigma <- exp(pars[3])
  L <- pars[4] ; U <- pars[5] ; n <- length(x)
  
  if (U<L) stop("Upper limit smaller than lower limit")
  
  # obtain the samples:
  mean <- x*gra + int
  sampled.size <- rnorm(n, mean, sigma)
  
  # perform rejection sampling to ensure the samples are within (L, U):
  while (any(sampled.size<L | sampled.size>U)){
    badSamples <- sampled.size<L | sampled.size>U
    replacements <- rnorm(sum(badSamples), mean[badSamples], sigma)
    sampled.size[badSamples] <- replacements
  }
  
  return(sampled.size)
}

sampleNorm <- function(x, pars){
  # purpose : produces samples as DTN without the double truncation.
  
  # extract parameter values:
  int <- pars[1] ; gra <- pars[2] ; sigma <- exp(pars[3])
  
  # obtain the samples:
  mean <- x*gra + int
  sampled.size <- rnorm(length(x), mean, sigma)
  return(sampled.size)
}

sampleOffNum <- function(z, rate){
  # purpose : Samples the number of offspring an individual has when they 
  #           reproduce, from a truncated poisson distribution
  # input   : The mean number of children each parent has when they reproduce
  # output  : The integer number of children had by the parent
  
  if (class(rate)!="numeric") stop("invalid rate input")
  rate <- exp(rate)
  
  # rejection algorithm that samples from the truncated Poisson with at least
  # 1 offspring produced:
  num <- rpois(length(z), rate)
  while(any(num==0)){
    badSamples <- num==0
    replacements <- rpois(sum(badSamples), rate)
    num[badSamples] <- replacements
  }
  
  return(num)
}

sampleExp <- function(z, pars){
  # purpose : Samples the size of an offspring from a truncated exponential
  #           distribution.
  # input   : The rate of the exponential distribution
  # output  : Samples of the children sizes
  
  rate <- exp(pars[1]) ; L <- pars[2] ; U <- pars[3]
  if (U<=L) stop('invalid size range')
  if (class(rate)!="numeric") stop("invalid rate input")
  if (rate<=0) stop("rate must be positive")
  
  # rejection algorithm that samples from the truncated exponential offspring
  # size distribution:
  samples <- rexp(length(z), rate)
  while(any(samples<L | samples>U)){
    badSamples <- samples < L | samples > U
    replacements <- rexp(sum(badSamples), rate)
    samples[badSamples] <- replacements
  }
  
  return(samples)
}


simulateIBM <- function(n, t, survFunc, survPars, growthSamp, growthPars,
                        reprFunc, reprPars, offNum=NULL, offNumSamp=NULL,
                        offNumPars=NULL, offSizeSamp, offSizePars, Start,
                        Schild=1, OneGend=FALSE, thresh=Inf,
                        CalcOffSurv=TRUE, popPrint=TRUE, verbose=FALSE){
  # purpose : Uses an Individual Based Model (IBM) to simulate census data used
  #           in fitting Integral Projection Models.
  # inputs  : n           - The starting size of the population
  #           t           - The maximum number of censuses to simulate
  #           survFunc    - The function which determines the probability of 
  #                         survival. It takes two arguments, the size of the 
  #                         individual, and a vector of parameters.
  #           survPar     - The vector of parameters for the survival function
  #           growthSamp  - The function which samples new sizes for
  #                         individuals. It takes two arguments, the size of the 
  #                         individual, and a vector of parameters
  #           growthPars  - The parameters for the function which samples the
  #                         new size of individuals.
  #           reprFunc    - The function which determines the probability of 
  #                         reproduction. It takes two arguments, the size of
  #                         the individual, and a vector of parameters.
  #           reprPars    - The parmeters for the reproduction function.
  #           offNum      - If not NULL, indiciates the fixed number of 
  #                         offspring that a parent has when they reproduce.
  #           offNumSamp  - If offNum is NULL, this is the function which
  #                         samples the number of offspring that a parent has
  #                         when they reproduce. It takes two arguments, the
  #                         size of the individual, and a vector of parameters
  #           offNumPars  - If offNum is NULL, this is the vector of parameters
  #                         for the number of offspring each parent has when 
  #                         they reproduce.
  #           offSizeSamp - The function which samples the size of the offspring
  #                         when they are first censused. It takes two
  #                         arguments, the size of the parent, and a vector of
  #                         parameters.
  #           offSizePars - The vector of parameters for the function which
  #                         samples the sizes of offspring.
  #           Start       - The size of the individual which fathers all the 
  #                         initial members of the population. Can be a vector,
  #                         if so, it must be of length n
  #           Schild      - The probability of survival of children to their 
  #                         first census.
  #           OneGend     - If TRUE, the census only tracks one gender, and so 
  #                         only half of the children born are assimilated into
  #                         the census data.
  #           thresh      - If the size of the population ever exceeds this 
  #                         quantity, we stop simulating.
  #           CalcOffSurv - If TRUE will calculate the number of offspring
  #                         that recruited for each parent in the population.
  #                         if FALSE, the function will fill this column with 
  #                         all NAs to save on computation time.
  #           popPrint    - If TRUE, the size of the population is printed at
  #                         the end of each census.
  #           verbose     - if TRUE, will print the survival rate, reproduction
  #                         rate, avergae change in size, rate of births and
  #                         child cenus rate at every survey:
  
  # make sure the number of children per litter is specified correctly:
  if (is.null(offNum) & (is.null(offNumSamp) | is.null(offNumPars))){
    stop('offNum or both offNumFunc and offNumPars must be specified')
  }
  
  if (is.null(thresh)) thresh <- .Machine$integer.max
  
  # to allow names ot be passed:
  survFunc %<>% match.fun ; growthSamp %<>% match.fun ; reprFunc %<>% match.fun
  offSizeSamp %<>% match.fun
  
  # If a parent always has the same number of children, we make sure it is 
  # an integer:
  if (!is.null(offNum)) offNum <- ceiling(offNum)
  
  # generate first generation from parent of fixed size:
  if (length(Start)!=n) Start <- rep(Start[1], n)
  initial.sizes <- offSizeSamp(Start, offSizePars)
  
  # generate the data.frame to store simulated data:
  DF <- data.frame(individual=1:n, size=initial.sizes, survived=1, 
                   census.number=1, parent.size=Start, reproduced=NA, 
                   prev.size=NA, off.survived=NA, off.born=NA)
  
  # create a data.frame which contains only the individuals from the previous
  # census (and will be constantly updated):
  previous.census <- DF
  
  time <- 2 # set the current census number (will be continually updated)
  while((time < t + 1) & nrow(previous.census)<thresh){
    
    # Select only the survivors of the previous census:
    survivorsDF <- subset(previous.census, previous.census$survived==1)
    current.pop.size <- nrow(survivorsDF)
    
    # Determine who survives this time:
    survProbs <- survFunc(survivorsDF$size, survPars)
    survivors <- rbinom(current.pop.size, 1, survProbs)
    if (sum(survivors)==0) break
    
    # Determine their new sizes:
    newSizes <- growthSamp(survivorsDF$size, growthPars)
    newSizes[-which(survivors==1)] <- NA # overwrite for dead individuals
    
    # Determine who reproduces (out of the survivors):
    reproducers <- rep(NA, current.pop.size)
    reprProbs <- reprFunc(na.omit(newSizes), reprPars)
    reproducers[survivors==1] <- rbinom(sum(survivors==1), 1, reprProbs)
    nRepr <- reproducers %>% na.omit %>% `==`(1) %>% sum
    parentsDF <- subset(survivorsDF, reproducers==1 & !is.na(reproducers))
    
    # Determine the number of offspring for each parent:
    if (!is.null(offNum)) censusOffNum <- rep(offNum, nRepr)
    else censusOffNum <- offNumSamp(parentsDF$size, offNumPars)
    
    # Determine the sizes of the offspring:
    reprSizes <- survivorsDF$size[reproducers==1 & !is.na(reproducers)]
    parentSizes <- rep(reprSizes, censusOffNum)
    parentIDs <- rep(parentsDF$individual, censusOffNum)
    offspringSizes <- offSizeSamp(parentSizes, offSizePars)
    
    # Determine the survival of the offspring:
    nOffspring <- length(offspringSizes)
    survivingChildren <- rbinom(nOffspring, 1, Schild)
    
    # Determine the sex of the offspring:
    probGend <- ifelse(isTRUE(OneGend), 1/2, 1)
    gender <- rbinom(nOffspring, 1, probGend)
    
    # Find out which children make it to census:
    newIDStart <- (max(DF$individual)+1)
    censusedChildren <- which(gender==1 & survivingChildren==1)

    # Print out summaries of the simulation if desired:
    if (verbose){
      cat("survival rate is", sum(survivors)/current.pop.size,"\n")
      oldSizes <- survivorsDF$size[which(survivors==1)]
      cat("average growth is", mean(na.omit(newSizes) - oldSizes), "\n")
      cat("reproduction rate is", nRepr/current.pop.size, "\n")
      cat("rate of births is", nOffspring/current.pop.size, "\n")
      cat("child censusing rate is", length(censusedChildren)/nOffspring,"\n")
      gr <- (length(censusedChildren)+sum(survivors))/current.pop.size
      cat("growth rate is", gr, "\n")
      cat("\n")
    }
    
    # Calculate how many surviving children per parent with a helper:
    helper <- function(x){
      ifelse(!x %in% parentIDs,
             NA,
             sum(x==parentIDs & survivingChildren==1))
    }#helper
    
    helper2 <- function(x){
      ifelse(!x %in% parentIDs,
             NA,
             sum(x==parentIDs))
    }#helper2
    
    # Only use the helper if the user wants us to. Otherwise make a vector of
    # NAs for speed of calculation:
    if (isTRUE(CalcOffSurv)){
      off.survived <- sapply(survivorsDF$individual, helper)
      off.born <- sapply(survivorsDF$individual, helper2)
    }
    
    else{
      off.survived <- rep(NA, current.pop.size)
      off.born <- rep(NA, current.pop.size)
    }
    
    # Create the DF for the parents:
    currentDF <- data.frame(individual=survivorsDF$individual,
                            size=newSizes,
                            survived=survivors,
                            census.number=time,
                            parent.size=NA,
                            reproduced=reproducers,
                            prev.size=survivorsDF$size,
                            off.survived=off.survived,
                            off.born=off.born)
    
    # Update the previous.census data.frame:
    if (length(censusedChildren)==0) previous.census <- currentDF
    
    else{
      newIDs <- newIDStart:(length(censusedChildren)+newIDStart-1)
      offspringDF <- data.frame(individual = newIDs,
                                size = offspringSizes[censusedChildren],
                                survived = 1, census.number = time,
                                parent.size = parentSizes[censusedChildren],
                                reproduced = NA, prev.size = NA,
                                off.survived = NA, off.born = NA)
      
      # Update the previous.census data.frame:
      previous.census <- rbind(currentDF, offspringDF)
    }
    
    # Update the overall data:
    DF <- rbind(DF, previous.census)
    
    # Iterate the census.number:
    time <- time + 1
    
    # Let the user know the size of the population:
    if(popPrint) print(nrow(previous.census))
  }
  
  return(DF)
}

# UNCOMMENT the below sections for example datasets to be generated:

# Try to simulate with similar values as the MSc Thesis from Duke:
# set.seed(102938)
# testSim <- simulateIBM(n=50, t=100,
#                        # set survival details:
#                        survFunc = linLogit, survPars = c(2.26, 0.23),
#                        # set growth details:
#                        growthSamp = sampleDTN,
#                        growthPars = c(0, 1, log(0.5), 0, 10),
#                        # set reproduction details:
#                        reprFunc = linLogit, reprPars = c(-3.5, 0.25),
#                        # set offspring number and size distribution details:
#                        offNumSamp = sampleOffNum, offNumPars = log(4),
#                        offSizeSamp = sampleExp,
#                        offSizePars = c(log(10), 0, 10),
#                        Schild=0.873,
#                        # set other miscelaneous parameters:
#                        Start=3, thresh=1000, OneGend = TRUE, popPrint = FALSE)
# 
# max.cens <- testSim %>% `$`(census.number) %>% max
# subset(testSim,  testSim$census.number==max.cens) %>% nrow %>% print

# # Look at the last census to give an idea on how the population grows:
# statement <- testSim$census.number==tsteps & testSim$survived==1
# statement %>% subset(x=testSim) %>% nrow %>% print

# Try to simulate Sheep data from 2.6 of Ellner, Childs & Rees:
# simmedData <- simulateIBM(n=500, t=1000,
#                           # set survival details:
#                           survFunc = linLogit, survPars = c(-9.65, 3.77),
#                           # set growth details:
#                           growthSamp = sampleDTN,
#                           growthPars = c(1.41, 0.56, log(0.08), 1.5, 3.55),
#                           # set reproduction details:
#                           reprFunc = linLogit, reprPars = c(-7.23, 2.6),
#                           # set offspring number and size distribution details:
#                           offNum=1, offSizeSamp = sampleDTN,
#                           offSizePars = c(0.36, 0.71, log(0.16), 1.5, 3.55),
#                           # Child survival probability:
#                           Schild=0.873,
#                           # set other miscelaneous parameters:
#                           Start=3, thresh=5000, OneGend = TRUE, popPrint = T)
