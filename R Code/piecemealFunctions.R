######################## PIECEMEAL SOLVING FUNCTIONS ###########################

# It'll be useful for us to have access to the size of the individual during the
# previous census for each observation, since this is the size which is relevant
# when attempting to estimate surival probability, reproduction probability etc:

getPreviousSizes <- function(DF){
  # purpose : gets the previous size of each individual in the data.frame 
  #           (with the previous size being the size of the individual at the 
  #           previous census)
  # inputs  : DF - The data.frame containing the IBM data
  # output  : A vector of length equal to the number of rows of DF, containing
  #           the previous size of each observation.
  
  N <- nrow(DF)
  prev.size <- rep(NA, N)
  
  for (i in 1:N){
    ID <- DF$individual[i]
    censN <- DF$census.number[i]
    if (censN==1) next # The first census won't have any 'previous size'
    previous.obs <- subset(DF, DF$individual==ID & DF$census.number==(censN-1))
    if (previous.obs$size %>% length %>% `==`(0)) next # skip if offspring
    prev.size[i] <- previous.obs$size
  }
  return(prev.size)
}

DesignMatrix = function(DataFrameInput, formula){
  # purpose : Produces the design matrix given a specified covariate formula, 
  #           and the data frame which contains all the required columns
  # inputs  : DataFrameInput - The data.frame which contains all of the 
  #                            columns which are specified in 'formula'
  #           formula        - The R formula object which specifies the 
  #                            relationship between the response variable and 
  #                            the covariates, all of which must be in 
  #                            'DataFrameInput'.
  # output  : A matrix object
  
  # We must remove the LHS of the formula, if it exists (it usually does):
  RHS = try(formula[[3]], silent=T)
  if (class(RHS)!='try-error'){ # This means there is a LHS as well as a RHS
    formula[[2]] = NULL         # This deletes the LHS and fixes the indices
  }
  
  # Make the design matrix for the whole data set:
  DM = model.matrix(object = formula, DataFrameInput)
  return(DM)
}

LinPredictor = function(parameters, DM){
  # purpose : Calculates the linear predictor for the parameter of interest,
  #           given the parameters and the design matrix
  # inputs  : parameters - The parameters for the model
  #           DM         - The design matrix for the model, as produced by
  #                        the 'DesignMatrix' function.
  # output  : A vector of real values, equal to the value of the linear 
  #           predictor for each observation.
  
  # Multiply DM by the column matrix of covariate parameters to obtain 
  # the linear predictor for the theta parameter of interest:
  parameters = matrix(parameters, ncol=1)
  
  # Exception handle the matrix multiplication to ensure the correct number 
  # of parameters has been supplied for the relationship specified by the 
  # passed formula object:
  LP = try({DM%*%parameters},silent=T)
  if (class(LP)=='try-error'){
    stop('Incorrect number of parameters supplied to LinPredictor')
  }
  return(LP)
}

reproductionNLL <- function(par, DF, formula=offspring~size){
  # purpose : Evaluates the negative log likelihood of the reproduction model, 
  #           when this model determines a probability of reproduction, and 
  #           a truncated Poisson distribution number of children when the 
  #           individual does reproduce.
  # inputs  : par - The parameters of the model. Should include the probability
  #                 of reproducing as the first entry, and the mean of the
  #                 truncated poisson distribution of the number of offspring
  #                 as the second entry. In the case where covariates are used, 
  #                 the 3:length(par) parameters should all be the covariates, 
  #                 and the first entry should be the intercept for the 
  #                 probability of reproducing (the covariates affect the
  #                 the probability of reproducing, not the mean number of
  #                 offspring).
  #          DF  - The data frame containing the data of the model. Should
  #                contain a column 'prev.size' for the data, as created by the
  #                getPreviousSizes function.
  # 
  # output  : A real number, the negative log likelihood of all the observations
  
  # remove the NA observations for offspring, since these represent just born
  # children, which can't have children yet (NOTE: will need to remove all 
  # individuals which are younger than sexual maturity when the model becomes
  # more complicated):
  DF <- subset(DF, !is.na(DF$offspring))
  x <- DF$offspring %>% as.numeric
  
  # Calculate the linear predictor for the probability of producing offspring:
  DM <- DesignMatrix(DF, formula)
  p <- LinPredictor(par[-2], DM) %>% plogis
  
  # log link for positive rate:
  rate <- exp(par[2])
  
  # The likelihood is split into two cases (by the law of total probability):
  # The case where the individual reproduces, and the case where it doesn't.
  firstTerm <- ((x>0)*p*dpois(x, rate))/(1-ppois(0, rate))
  secondTerm <- (x==0)*(1-p)
  (firstTerm + secondTerm) %>% log %>% sum %>% `-` %>% return
}

growthNLL <- function(par, DF, L, U, formula=~prev.size){
  # purpose : Evluates the negative log likelihood of the growth function, 
  #           which is a doubly-truncated normal distribution with mean
  #           typically centred at a function of the size of the individual at
  #           the previous census.
  # inputs  : par     - The vector of parameters. The first parameter should be
  #                     the intercept of the mean parameter, the second should
  #                     be the standard deviation of the truncated normal 
  #                     distribution, and all remaining parameters should be 
  #                     the required covariates specified in the formula for the
  #                     expectation of the size of the individual at the current
  #                     census
  #           DF      - The data.frame containing all of the required variables
  #                     specified by 'formula', and must contain a 'size' 
  #                     column with the current size of the individual 
  #           L       - The lower limit of the size variable
  #           U       - The upper limit of the size variable
  #           formula - The formula which specifies the relationship which
  #                     should be used to calculate the expected value of the
  #                     size of the individual at the current census.
  # output  : A real number, the negative log likelihood of all the observations
  # DF <- subset(DF, !is.na(DF$size) & DF$census.number>1)
  if (any(DF$size<L | DF$size>U)) stop("contains invalid size values")
  
  # extract the response variable using the formula:
  isLHS <- try(formula[3], silent=T)
  if (class(isLHS)=='try-error') stop('invalid formula, LHS must be specified')
  LHS <- formula[2] %>% as.character %>% sub(pattern="()",replacement="")
  response <- DF[LHS][[1]] %>% as.numeric
  
  growth.pars <- par[-3]
  sigma <- par[3] %>% exp # log link for positive SD
  
  # Make the linear predictor:
  DM <- DesignMatrix(DF, formula)
  mean <- LinPredictor(growth.pars, DM)
  
  # return the likelihood:
  denom <- log(pnorm(U, mean, sigma) - pnorm(L, mean, sigma))
  dnorm(response, mean, sigma) %>% log %>% `-`(denom) %>% sum %>% `-` %>% return
}

growthLinear <- function(zPrime, z, par, L, U, lik=FALSE){
  # purpose : A form of the growth function with a linear relationship between
  #           previous size and size at the current census
  # inputs  : zPrime - The size at the current census
  #           z      - The size at the previous census
  #           par    - The parameters of the doubly truncated normal
  #                    distribution of size given previous size
  #           L      - The lower limit of the size variable
  #           U      - The upper limit of the size variable
  #           lik    - If TRUE, returns the product of all of the observations,
  #                    so that the output is the likelihood
  # output  : A single number, the liklihood of the observed size, given the
  #           previous size
  mean <- par[1]+ par[2]*z
  sigma <- par[3] %>% exp
  denom <- pnorm(U, mean, sigma) - pnorm(L, mean, sigma)
  calcs <- dnorm(zPrime, mean, sigma)/denom
  if (isTRUE(lik)) calcs <- prod(calcs)
  return(calcs)
}

# A quick function for when a model holds a quantity constant:
returnConstant <- function(x, const) const

offSizeNLL <- function(par, L, U, DF){
  # purpose : Calculates the negative log likelihood of the exponential 
  #           recruitment function for all newborns in the dataset.
  # inputs  : par - The rate parameter for the exponential distribution
  #           L   - The lower limit for the size variable
  #           U   - The upper limit for the size variable
  #           DF  - The data.frame containing the data for the analysis
  # output  : A single real number, the value of the negative log likelihood
  #           evaluated for the sizes of all offspring
  
  # scan through the data to select only the children on the census that they
  # were born on:
  DF <- subset(DF, is.na(DF$offspring) & DF$survived==1 & is.na(DF$prev.size))
  y <- DF$size
  
  if (any(y<L | y>U)) stop('invalid offspring sizes')
  rate <- exp(par)
  numer <- dexp(y, rate)
  denom <- pexp(U, rate) - pexp(L, rate)
  (numer/denom) %>% log %>% sum %>% `-` %>% return
}

recruitmentLik <- function(z, par, L, U){
  # purpose : calculates the likelihood of an exponential recruitment function
  #           given a vector of sizes rather than an input data.frame
  # inputs  : z   - The scalar or vector sizes of the individuals when they 
  #                 are first censused
  #           par - The rate parameter of the exponential distribution
  #           L   - The lower limit of the size variable
  #           U   - The upper limit of the size variable
  if (any(z<L | z>U)) stop('invalid offspring sizes')
  rate <- exp(par)
  numer <- dexp(z, rate)
  denom <- pexp(U, rate) - pexp(L, rate)
  (numer/denom) %>% prod %>% return
}

doublyTruncatedNormal <- function(y, x, pars){
  # purpose : Calculates the probability of an individual having size y given 
  #           it was size x at the previous time step. Has the form of a 
  #           doubly truncated normal with a mean parameter
  # inputs  : y         - The size of the individual at the current time step 
  #           x         - The size of the individual at the previous time step
  #           pars:
  #             intercept - The intercept for the mean parameter
  #             gradient  - The gradient for the mean parameter
  #             sigma     - The standard deviation of the normal distn for the
  #                         size of the individual at the current time step, 
  #                         given the size at the previous time step
  #             L         - The lower bound of the size variable
  #             U         - The upper bound of the size variable
  # output : The growth density evaluated at x and y
  
  # Extract parameters from the passed vector:
  intercept <- pars[1] ; gradient <- pars[2] ; sigma <- pars[3] %>% exp
  L <- pars[4] ; U <- pars[5]
  
  if (!L<U) stop("Growth: Upper bound does not exceed Lower bound")
  if (any(x<L) | any(x>U) | any(y<L) | any(y>U)) stop("Growth: Invalid observation")
  
  mean <- intercept + x*gradient          # Calculate the mean
  p <- dnorm(y, mean, sigma)              # Calculate the probabilty of y|x
  lowerTail <- pnorm(L, mean, sigma)      # Lower tail prob
  upperTail <- pnorm(U, mean, sigma)      # Upper tail prob
  denom <- upperTail - lowerTail          # Truncate the dist on both tails
  return(p/denom)
}

normal <- function(y, x, pars){
  # Extract parameters from the passed vector:
  intercept <- pars[1] ; gradient <- pars[2] ; sigma <- pars[3] %>% exp
  
  mean <- intercept + x*gradient          # Calculate the mean
  dnorm(y, mean, sigma)  %>% return       # Calculate the probabilty of y|x
}

kernelOneVar <- function(m, growthFunc, growthPars, survFunc, survPars,
                         repFunc, repPars, offNum, offSizeFunc, offSizePars, L, 
                         U, childSurv=1, halfPop=F){
  # purpose : Evaluates the Integral Projection Model Kernal for one time step,
  #           given the range of the size values, and the number of points to 
  #           use with the midpoint.
  # inputs  : m          - The dimension of the iteration matrix
  #           growthFunc - The R function which evaluates the probability 
  #                        of growing from z to zPrime
  #           growthPars - The parameters of the growth function
  #           survFunc   - The R function which determines the probability of 
  #                        survival of the individual
  #           survPars   - The parameters of the survival function
  #           repFunc    - The function which determines the probability that
  #                        an individual has of reproducing
  #           repPars    - The parameters of the reproduction function
  #           offNum     - The mean number of offspring produced by each 
  #                        reproducing individual
  #           recrFunc   - The R function which determines the distribution of
  #                        of sizes of newborn individuals
  #           recrPars   - The parameters of the recruitment function
  #           childSurv  - The probability each child has of surviving till the
  #                        next census
  #           halfPop    - If TRUE, this indicates that the data set only
  #                        keeps track of one sex, which is assumed to make up
  #                        roughly half of the population
  #           offSizeZp  - If TRUE, the size of the offspring depends on the
  #                        current size of the individual. The size of the
  #                        offspring is assumed to be the first parameter.
  # output : A single number, the evaluated Kernel.
  #
  # Note: uses formula 2.3.6 of Ellner, Childs & Rees for the kernel of a one
  #       variable post-reproductive IPM.
  
  # produce mesh points for the midpoint rule:
  h <- (U-L)/m
  meshpts <- L + (1:m)*h - h/2
  
  # To account for only tracking half the population:
  multiplier <- ifelse(is.null(halfPop), 1, 1/2)
  
  # calculate first term of the kernel:
  pFunc <- function(zprime, z){
    survivalProb <- survFunc(z, survPars)
    growthProb <- growthFunc(zprime, z, growthPars)
    return(survivalProb*growthProb)
  }
  
  # Calculate the second term of the kernel:
  fFunc <- function(zprime, z){
    survivalProb <- survFunc(z, survPars)
    repProb <- repFunc(z, repPars)
    offSizeProb <- offSizeFunc(zprime, z, offSizePars)
    return(survivalProb*repProb*offNum*offSizeProb*multiplier*childSurv)
  }
  
  Pterm <- h*outer(meshpts, meshpts, pFunc)
  Fterm <- h*outer(meshpts, meshpts, fFunc)
  kernel <- Pterm + Fterm
  rownames(kernel) <- meshpts
  colnames(kernel) <- meshpts
  class(kernel) <- "IterationMatrix"
  return(kernel)
} 
