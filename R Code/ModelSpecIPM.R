# THIS FILE SPECIFIES THE PRIOR FUNCTION AND THE LOG TARGET FOR THE POSTERIOR ##
# OF THE 'PARTICLE FILTER' VERSION OF THE IPM. #################################

evalPriors <-  function(vals, funcs, listOfPars, Log=T, checks=F){
  # purpose : Evaluates the log likelihood of all the priors.
  # inputs  : vals       - The vector of sampled parameter values.
  #           funcs      - The vector of character names for the functions which
  #                        evaluate the priors.
  #           listOfPars - The list of same length as vals and funcs, which 
  #                        contains other arguments which should be passed to
  #                        functions which evaluate the priors. This is a list
  #                        of lists, as each parameter may have a prior with 
  #                        multiple parameters.
  #           log        - If TRUE, calculates the log of the prior. If FALSE,
  #                        it is assumed that these priors already return a log
  #                        probability.
  #           checks     - if TRUE, will complete some basic checks on the
  #                        arguments of this function.
  
  if (checks){
    if (class(vals)!='numeric') stop('Invalid vals input')
    if (class(funcs)!='character') stop('Invalid funcs input')
    if (class(listOfPars)!='list') stop('Invalid listOfPars input')
    if (class(log)!='logical') stop('Invalid log input')
    n <-  length(vals) ; y <- length(funcs) ; z <- length(listOfPars)
    if (n!=y | n!=z | y!=z) stop('Incompatible argument dimensions')
  }
  
  # Avoids a bug:
  names(vals) <- NULL
  
  # Define the helper function which gets the answer we require:
  h <- function(index) do.call(match.fun(funcs[index]),
                               c(vals[index], listOfPars[[index]]))
  
  # setup the log operation if required:
  if (!Log) logfunc <- function(x) x
  else logfunc <- log
  
  # Loop through the priors for all parameters
  sapply(1:length(vals), h) %>% logfunc %>% sum %>% return
}

logTargetIPM <- function(proposed, logTargetPars, returnNeg = F, check = F, 
                         returnW = F, printProp = F, returnLL = F){
  # purpose : Evaluates the log target distribution of the IPM state space model
  #           given a list of parameters which defines the specification: 
  # inputs  : proposed      - A vector of new proposed parameter values
  #           returnNeg     - If TRUE, returns the negative log target instead
  #           check         - If TRUE, does a check for some inputs
  #           returnW       - If TRUE, will return the standardised weights of
  #                           the particle filter instead
  #           returnLL      - If TRUE, returns the Log likelihood instead of
  #                           the log posterior
  #           logTargetPars - The required parameters for this specific model
  #                           this parameter is a list which contains:
  #                           
  #                           priorFunc - A function which evaluates the prior
  #                                       probability of a proposed parameter
  #                                       set. It is assumed that this function
  #                                       takes as its first input the parameter
  #                                       set, and all inputs specified by
  #                                       priorPars afterwards.
  #                           priorPars - The parameters for each of the prior
  #                                       distributions.
  #                           skeleton  - A list which allows us to extract the
  #                                       parameters in proposed, without 
  #                                       imposing a fixed structure on the 
  #                                       required number of parameters per 
  #                                       vital rate function.
  #                           survFunc  - The function which evaluates the 
  #                                       probability of survival
  #                           growthSamp - The function which samples the growth
  #                                        of individuals.
  #                           reprFunc  - The function which determines the 
  #                                       probability of survival.
  #                           offNumSamp - The function which samples the
  #                                        number of offspring a parent has.
  #                           offSizeSamp - The function which samples the 
  #                                         sizes of newborns.
  #                           oneSex    - Optionally, specifies if the survey
  #                                       is only tracking one sex of a
  #                                       population which is assumed to birth
  #                                       equal proportions of each sex. Default
  #                                       value is TRUE.
  #                           shift     - Optionally, shift can be passed as a
  #                                       fixed value rather than a parameter.
  #                                       see sampleStateIPM for information on
  #                                       'shift'.
  #                           mu        - The function which allows us to draw
  #                                       samples from the initial distribution 
  #                                       of size class counts.
  #                           muPar     - The parameters of the initial
  #                                       distribution of size class counts.
  #                           b         - The number of samples used by the 
  #                                       particle filter to estimate the log
  #                                       likelihood.
  #                           Y         - The matrix of count observations. Each
  #                                       column is a time step.
  #                           obsProb   - The function which evaluates the
  #                                       likelihood of the observations given 
  #                                       the hidden state
  #                           breaks    - The vector of breakpoints which 
  #                                       gives the range of each size class
  #
  # output  : A single real number, the log posterior for the given parameter
  #           values.
  #
  # NOTE : It is the user's responsibility to double check that the parameters
  #        being fitted, and their order in the skeleton is equal to the exact
  #        order in which they have specified the priors.
  if (check){
    if(length(logTargetPars$muPar[[2]])!=(length(logTargetPars$breaks)-1)){
      stop('Dimension mismatch in initial distribution and size class numbers')
    }
  }
  
  if (returnNeg) multiplier <- -1
  else multiplier <- 1
  if (printProp) print(proposed)
  # Get the parameters into a list:
  skeleton <- logTargetPars$skeleton
  vectorPars <- proposed
  proposed %<>% relist(skeleton=skeleton)
  
  # necessary checks:
  if (is.null(proposed$shift) & is.null(logTargetPars$shift)) stop('need shift')
  
  # Extract the function parameters:
  survPars <- proposed$survPars
  growthPars <- proposed$growthPars
  reprPars <- proposed$reprPars
  offNumPars <- proposed$offNumPars
  offSizePars <- proposed$offSizePars
  Schild <- proposed$Schild
  shift <- proposed$shift ; if (is.null(shift)) shift <- logTargetPars$shift
  obsProbPar <- proposed$obsProbPar
  
  # Extract function arguments:
  b <- logTargetPars$b
  Y <- logTargetPars$Y
  obsProb <- logTargetPars$obsProb %>% match.fun
  survFunc <- logTargetPars$survFunc
  growthSamp <- logTargetPars$growthSamp
  reprFunc <- logTargetPars$reprFunc
  offNumSamp <- logTargetPars$offNumSamp
  oneSex <- logTargetPars$oneSex
  if (is.null(logTargetPars$oneSex)) oneSex=T
  offSizeSamp <- logTargetPars$offSizeSamp
  
  # Extract prior arguments:
  mu <- logTargetPars$mu
  muPar <- logTargetPars$muPar
  priorFunc <- logTargetPars$priorFunc %>% match.fun
  priorPars <- logTargetPars$priorPars
  breaks <- logTargetPars$breaks
  funcs <- logTargetPars$priors
  
  # Create the list of arguments required for the state space sampler:
  stateSpaceSampArgs <- list(survFunc = survFunc, survPars = survPars,
                             growthSamp = growthSamp, growthPars = growthPars,
                             reprFunc = reprFunc, reprPars = reprPars, 
                             offNumSamp = offNumSamp, offNumPars = offNumPars,
                             offSizeSamp = offSizeSamp, breaks = breaks,
                             offSizePars = offSizePars, Schild=Schild,
                             oneSex = oneSex, shift = shift)
  
  # Get an estimate of the log likelihood from the particle filter:
  ll <- particleFilter(Y=Y, mu=mu, muPar=muPar, obsProb = obsProb,
                       sampleState = vectorisedSamplerIPM,
                       sampleStatePar = stateSpaceSampArgs,
                       obsProbPar = obsProbPar, b = b, returnW = returnW)
  
  if (returnW) return(ll)
  
  else{
    if (returnLL) return(ll)
    else # Get the evaluation of the priors:
      (ll + priorFunc(vectorPars, funcs, priorPars)) %>%
        `*`(multiplier) %>% return
  }
}

posteriorGrowthRate <- function(chains, IPMLTP, growthFunc, offSizeFunc, L=0, U,
                                level = 0.95, m = 500, printRate = NULL,
                                cluster, clNeeds = NULL){
  # purpose : produces a credible interval for the growth rate of the population
  #           using the posterior samples from an MCMC chain.
  # inputs  : chains      - The list of matrices which contain the samples from
  #                         the posterior produced by the MCMC chains.
  #           IPMLTP      - The logTargetParameters from the call to the pMH
  #                         function, which specifies the IPM which was used.
  #           growthFunc  - The function which evaluates the probability of
  #                         growth.
  #           offsizeFunc - The function which evaluates the probability of a 
  #                         child's size given the size of its parent.
  #           L           - The lower limit of the size variable.
  #           U           - The upper limit of the size variable.
  #           level       - The percentage of the probability mass of the
  #                         posterior distribution of growth rate that the
  #                         credible interval should cover.
  #           m           - The number of gridpoints which should be used to 
  #                         estimate the growth rate (the gridsize of the
  #                         numeric integration of the kernel).
  #           printRate   - Will print the number of the sample the function is
  #                         working through every printRate samples.
  #           cluster     - A cluster to run the function as a foreach loop in 
  #                         parallel
  growthFunc %<>% match.fun ; offSizeFunc %<>% match.fun
  quantiles <- c(1-level, 1+level)/2
  
  # extract the required information from IPMLTP:
  if (is.null(IPMLTP$oneSex)) oneSex=T
  else oneSex <- IPMLTP$oneSex
  skeleton <- IPMLTP$skeleton
  survFunc <- IPMLTP$survFunc %>% match.fun
  reprFunc <- IPMLTP$reprFunc %>% match.fun
  
  # combine the chains together:
  samples <- chains[[1]]
  for (i in 2:length(chains)){samples %<>% rbind(chains[[i]])}
  
  registerDoParallel(cluster)
  clusterExport(cluster, c("kernelOneVar", clNeeds))
  
  ls <- foreach(i = 1:nrow(samples)) %dopar% {
    if (!is.null(printRate) & i%%printRate==0) print(i)
    # turn the sample into a list:
    sample <- samples[i, -1]
    sample %<>% relist(skeleton=skeleton)
    
    # extract the parameters:
    survPars <- sample$survPars
    growthPars <- sample$growthPars
    reprPars <- sample$reprPars
    offNumPars <- sample$offNumPars
    offSizePars <- sample$offSizePars
    Schild <- sample$Schild
    
    # get the growth rate:
    to.ls <- kernelOneVar(m = m, growthFunc = growthFunc,
                               growthPars = growthPars, survFunc = survFunc,
                               survPars = survPars, repFunc = reprFunc,
                               repPars = reprPars, offNum = offNumPars,
                               offSizeFunc = offSizeFunc,
                               offSizePars = offSizePars, L = L, U = U,
                               childSurv = Schild, halfPop = oneSex) %>%
      eigen %>% `$`(values) %>% `[`(1) %>% Re
  }
  
  results <- unlist(ls)
  CI <- quantile(results, probs = quantiles)
  list(CI = CI, results = results) %>% return
}
