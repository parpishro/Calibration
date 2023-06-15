#' Specify priors for each categories of parameters
#'
#' `setup_prior` uses user-given prior types and parameters to select and parametrize
#' the prior function and compute its mean to be used as initial value. This function
#' setup prior functions to return log of prior distribution.
#'
#' @param prior A String that specifies the type of the prior
#' @param p1    A double as first parameter of the prior distribution
#' @param p2    A double as second parameter of the prior distribution
#' @noRd
#' @return A function that given its input, x, computes the log of the chosen prior probability density function
setup_prior <- function(prior, p1, p2) {

  log_prior <- switch (prior,
                       uniform      = function(x) log(dunif(x,  min      = p1, max    = p2)),
                       gaussian     = function(x) log(dnorm(x,  mean     = p1, sd     = p2)),
                       gamma        = function(x) log(dgamma(x, shape    = p1, scale  = p2)),
                       beta         = function(x) log(dbeta(x,  shape1   = p1, shape2 = p2)),
                       lognormal    = function(x) log(dlnorm(x, meanlog  = p1, sdlog  = p2)),
                       logistic     = function(x) log(dlogis(x, location = p1, scale  = p2)),
                       betashift    = function(x) log(dbeta(x-1, shape1   = p1, shape2 = p2)),
                       exponential  = function(x) log(dexp(x,   rate     = p1)),
                       inversegamma = function(x) log(((p2^p1)/gamma(p1)) * (x^(-p1-1)) * exp(-p2/x)),
                       jefferys     = function(x) -0.5*log(x),
                       stop("at least one prior missing or invalid!"))

  return(log_prior)
}


