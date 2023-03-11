#' Specify priors for each categories of parameters
#'
#' `setup_prior` uses user-given prior types and parameters to select and parametrize
#' the prior function and compute its mean to be used as initial value. This function
#' setup prior functions to return log of prior distribution.
#'
#' @param prior A String that specifies the type of the prior
#' @param p1    A double as first parameter of the prior distribution
#' @param p2    A double as second parameter of the prior distribution
#'
#' @return A function that given its input, x, computes the log of the chosen prior probability density function
setup_prior <- function(prior, p1, p2) {

  log_prior <- switch (prior,
                       uniform      = function(x) dunif(x,  min      = p1, max    = p2, log = T),
                       gaussian     = function(x) dnorm(x,  mean     = p1, sd     = p2, log = T),
                       gamma        = function(x) dgamma(x, shape    = p1, scale  = p2, log = T),
                       beta         = function(x) dbeta(x,  shape1   = p1, shape2 = p2, log = T),
                       lognormal    = function(x) dlnorm(x, meanlog  = p1, sdlog  = p2, log = T),
                       logistic     = function(x) dlogis(x, location = p1, scale  = p2, log = T),
                       exponential  = function(x) dexp(x,   rate     = p1,              log = T),
                       inversegamma = function(x) p1*log(p2) - lgamma(p1) - (p1+1)*log(abs(x)) - p2/(abs(x)),
                       jefferys     = function(x) -0.5*log(x),
                       logbeta      = function(x) -log(8)- x/4 - 0.5*log(1-exp(-x/4)),
                       stop("at least one prior missing or invalid!"))

  return(log_prior)
}
