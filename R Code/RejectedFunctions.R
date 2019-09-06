# This file serves as a graveyard for some of the rejected functions throughout
# this project, when belter alternatives were found through profiling or
# changing the model specification.

# slower than using an rmultinom version, so discarded.
uniformMuREJECTED <- function(popSize, traitClasses){
  # purpose : returns a vector of individual counts which is uniformly 
  #           distributed across the dimensions of the state space
  # inputs  : popSize      - The integer total population size (to be spread
  #                          uniformly amongst the train classes) 
  #           traitClasses - The integer number of bins that the trait is split
  #                          into for the purpose of the analysis.
  samp <- sample(1:traitClasses, popSize, replace=T)
  sapply(1:traitClasses, function(x) sum(samp==x)) %>% return
}

# rejected for being slightly slower than the versionusing lots of pipes:
sampleStateIPMREJECTED <- function(previousState, survFunc, survPars,
                           growthSamp, growthPars, reprFunc, reprPars, 
                           offNumSamp, offNumPars, offSizeSamp, offSizePars,
                           Schild, L, U, oneSex=TRUE, breaks=NULL,
                           checks=FALSE, sizeHist=FALSE, unif=FALSE){
  # purpose : A barebones version of simulateIBM which takes as input a size 
  #           distribution at a given time, and given the parameters of the 
  #           vital rate functions, produces a sample of the size distribution
  #           of the individuals at the next time step (census).
  # inputs  : previousState - A vector of the number of individuals in each
  #                           size class. Can be any length.
  #           survFunc      - The function which evaluates the probability of 
  #                           survival. The first argument is the size of the
  #                           individual.
  #           survPars      - The parameters for the survival function.
  #           growthSamp    - The function which samples the new size of
  #                           individuals, given their current size. The first
  #                           argument is the size.
  #           growthPars    - The parameters of the growth function.
  #           reprFunc      - The function which calculates the probability of
  #                           survival. The first argument is the size.
  #           reprPars      - The parameters for the reproduction function.
  #           offNumSamp    - The function which samples from the distribution
  #                           of the number of offspring. 
  #           offNumPars    - The parameters for the offspring number sampling
  #                           distribution.
  #           offSizeSamp   - The function which samples from the distribution
  #                           of offspring sizes given parent sizes. Parent size
  #                           is the first argument.
  #           offSizePars   - The parameters of the child size distribution.
  #           Schild        - The probability that a child survives to the 
  #                           first census.
  #           L             - The lower limit of the size distribution.
  #           U             - The upper limit of the size distribution.
  #           oneSex        - If TRUE, the survey is tracking only one sex, and
  #                           so we simulate the sex of the children with a 
  #                           Bernoulli(0.5) distribution before including them
  #                           in the survey.
  #           breaks        - The breakpoints of the size classes. Should be a 
  #                           vector of length D+1, where D is the number of
  #                           size classes. The extreme entries should be L and
  #                           U, but these will be replaced by -Inf and Inf
  #                           when producing the new size distribution, since
  #                           some non-truncated samplers for the growth and 
  #                           offspring size distributions will produce values
  #                           slightly outside of [L, U]
  #           checks        - If TRUE, performs some input checks
  #           sizeHist      - If TRUE, produces a histogram of the sizes used
  #                           to project the state forwards
  #           unif          - If TRUE, disperses points in size classes
  #                           uniformly in that class, to avoid systematic 
  #                           bias
  # output  : A vector of counts of the same length as previousState
  # note : Could add functionality for customised breaks, maybe useful if 
  #        population is very clustered in a certain range.
  # note : All function arguments can be passed as the character name of the
  #        function instead
  require(plyr)
  
  # To insure function names work too:
  survFunc %<>% match.fun ; growthSamp %<>% match.fun ; reprFunc %<>% match.fun
  offNumSamp %<>% match.fun ; offSizeSamp %<>% match.fun
  
  # if half the children are of a sex not being tracked:
  if (isTRUE(oneSex)) Schild %<>% `/`(2)
  
  # Now the implementation:
  D <- length(previousState)
  
  # Some small setup for discretisiation:
  if (is.null(breaks)){
    stepsize <- (U-L)/D
    breaks <- L + (0:D)*stepsize
    sizes <- L + (1:D)*stepsize - stepsize/2
  }
  
  # Get midpoint of each user defined interval:
  else{
    sizes <- breaks[-(D+1)] + diff(breaks)/2
  }
  
  # Uniformly disperse the simulated points throughout the size class to avoid
  # any kind of systematic bias (overwriting the old sizes):
  if (unif){
    Lvec <- rep(breaks[-(D+1)], previousState)
    Uvec <- rep(breaks[-1], previousState)
    sizes <- runif(sum(previousState), Lvec, Uvec)
  }
  
  else sizes %<>% rep(previousState)
  if (sizeHist) hist(sizes)
  
  # Survival and growth:
  survProbs <- survFunc(sizes, survPars)
  survivingSizes <- weightedSelection(sizes, survProbs)
  grownSizes <- growthSamp(survivingSizes, growthPars)
  
  # Reproduction and survival of the offspring:
  reproductionProbs <- reprFunc(grownSizes, reprPars)
  reproducedSizes <- weightedSelection(grownSizes, reproductionProbs)
  childNums <- offNumSamp(reproducedSizes, offNumPars)
  reproducedSizes <- rep(reproducedSizes, childNums)
  reproducedSizes <- weightedSelection(reproducedSizes, Schild)
  childSizes <- offSizeSamp(reproducedSizes, offSizePars)
  
  # So breaks include counts for any sizes:
  breaks[c(1, D+1)] <- c(-Inf, Inf)
  
  # Count the indivs in each size class:
  output <- findInterval(c(childSizes, grownSizes), breaks) %>% count
  zeroCounts <- try(data.frame(x=setdiff(1:D, output$x), freq=0), silent=T)
  if (class(zeroCounts)=='try-error') return(output$freq)
  else zeroCounts %>% rbind(output) %>% `$`(freq) %>% return
}

# works fine, replaced by Rcpp code:
vectorToCountsREJECTED <- function(vector, breaks){
  # purpose : A function which takes a vector and some breakpoints, and returns
  #           the number of values in vector that reside in each of the
  #           intervals defined by the breaks.
  # inputs  : vector - The vector of values we wish to count up
  #           breaks - The breaks which define the interval boundaries
  # output  : A vector of counts of length length(breaks) - 1.
  require(plyr)
  D <- length(breaks) - 1
  output <- findInterval(vector, breaks) %>% count
  zeroCounts <- try(data.frame(x=setdiff(1:D, output$x), freq=0), silent=T)
  if (class(zeroCounts)=='try-error') return(output$freq)
  else zeroCounts %>% rbind(output) %>% `$`(freq) %>% return
}

# works fine, replaced by Rcpp code:
weightedSelectionREJECTED <- function(x, probs){
  # purpose : Takes an input vector and keeps each element with the probability
  #           determined by probs
  # inputs  : x     - The vector from which we sample values
  #           probs - The probability of keep each element
  # output  : A vector of length <= length(x) containing the elements of x
  #           that weren't removed
  keeping <- rbinom(length(x), 1, probs) %>% as.logical %>% which
  return(x[keeping])
}

# works fine, replaced by Rcpp code:
countsToProbs <- function(projectionOutput){
  # purpose : Produces a vector of probabilities of being in a given size class
  #           given the output of a function which projects a state over many
  #           time steps
  # input   : projectionOutput - The output of a function such as
  #                              projectStateSpace
  # output  : A vector of probabilities of the same length as the number of 
  #           size classes used in the simulated forward projection resulting
  #           in projectionOutput.
  n <- ncol(projectionOutput)
  if (is.null(n)) return(projectionOutput/sum(projectionOutput))
  else return(projectionOutput[,n]/sum(projectionOutput[,n]))
}
