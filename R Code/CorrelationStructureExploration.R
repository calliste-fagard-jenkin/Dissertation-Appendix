library(magrittr)
source('piecemealFunctions.R') # for the DesignMatrix and LinPredictor funcs

rCNS <- function(df, formula1, formula2, pars1, pars2, sigmaE){
  # purpose : Simulates samples from two normal distributions, where the 
  #           expected values of each normal are linked by a common term epsilon
  #           which is drawn from a normal distribution of mean 0 and sd sigmaE.
  # inputs  : df       - The data frame containing all of the dependent and
  #                      independent variables required by the formulae.
  #           formula1 - The formula which calculates the expected value of the
  #                      first normal distribution
  #           formula2 - The formula which calculates the expected value of the
  #                      second normal distribution
  #           pars1    - The parameters for the first normal distribution. The
  #                      final parameter should be the standard deviation, the 
  #                      previous ones should be the parameters required to 
  #                      calculate the expected value of this distribution given 
  #                      formula1.
  #           pars2    - The parameters for the second normal distribution. The
  #                      final parameter should be the standard deviation, the 
  #                      previous ones should be the parameters required to 
  #                      calculate the expected value of this distribution given 
  #                      formula2.
  #           sigmaE   - The standard deviation of the common term that the 
  #                      expected value of the two normal distributions share.
  # output  : A list with two components, the first is the samples from the
  #           first normal, the second the samples from the second.
  if (sigmaE<=0 | pars1[length(pars1)]<=0 | pars2[length(pars2)]<=0) stop(
    'Standard deviation parameters must be positive'
  )
  
  # split the parameter vectors for clarity:
  sigma1 <- pars1[length(pars1)] ; sigma2 <- pars2[length(pars2)]
  pars1 <- pars1[-length(pars1)] ; pars2 <- pars2[-length(pars2)]
  
  # generate the deviates in a few steps:
  n <- nrow(df)
  epsilon <- rnorm(n, 0, sigmaE)
  mu1 <- DesignMatrix(df, formula1) %>% LinPredictor(parameters = pars1)
  mu2 <- DesignMatrix(df, formula2) %>% LinPredictor(parameters = pars2)
  deviates1 <- rnorm(n, mu1 + epsilon, sigma1)#*(sigma1/sqrt(sigma1^2+sigmaE^2))
  deviates2 <- rnorm(n, mu2 + epsilon, sigma2)#*(sigma2/sqrt(sigma2^2+sigmaE^2))
  return(list(X1=deviates1, X2=deviates2))
}

# By playing around changing the value of sigmaE we can see that the
# correlation between X1 and X2 increases, introducing the correlation between 
# the parameters that we want in the hierarchical bayes model. However this also
# increases the dispersion parameter, and so the parameters for the variances
# of X1 and X2 will be heavily correlated with sigmaE.
npoints <- 10000
df <- data.frame(size = rnorm(npoints, 3, 0.6), age = rpois(npoints, 50))
rCNS(df, ~size, ~age, c(0.5, 3, 0.4), c(-10, 1, 0.2), 20) %>% as.data.frame %>%
  plot(col=adjustcolor('black', 0.2), pch=16)
