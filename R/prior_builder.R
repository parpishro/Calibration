#' Building Prior Function
#'
#' `prior_builder` uses user-given prior types and parameters to select and parametrize
#' the prior function. This function
#' setup prior functions to return log of prior distribution.
#'
#' @param prior character string representing prior distribution family
#' @param p1    double representing the first parameter of the prior distribution
#' @param p2    double representing the second parameter of the prior distribution
#' @param ln    logical denoting whether the output prior function return log of value.
#'              Default is True as it log priors are needed in computations and they
#'              simplifies the computation.
#'
#' @return function that given its input, x, computes the log of the chosen prior
#'         probability density function
#' @export
prior_builder <- function(prior, p1, p2) {
  force(p1)
  force(p2)
  prior_fn <- switch (prior,
                      uniform      = function(x) return(dunif(x,   min=p1,      max=p2,    log=T)),
                      gaussian     = function(x) return(dnorm(x,   mean=p1,     sd=p2,     log=T)),
                      gamma        = function(x) return(dgamma(x,  shape=p1,    scale=p2,  log=T)),
                      beta         = function(x) return(dbeta(x,   shape1=p1,   shape2=p2, log=T)),
                      lognormal    = function(x) return(dlnorm(x,  meanlog =p1, sdlog=p2,  log=T)),
                      logistic     = function(x) return(dlogis(x,  location=p1, scale=p2,  log=T)),
                      betashift    = function(x) return(dbeta(x-1, shape1=p1,   shape2=p2, log=T)),
                      exponential  = function(x) return(dexp(x,    rate=p1,                log=T)),
                      inversegamma = function(x) return(p1*log(p2) - log(gamma(p1)) - (p1+1)*log(x) - p2/x),
                      jefferys     = function(x) return(-0.5*log(x)),
                      stop("at least one prior missing or invalid!"))
  return(prior_fn)
}


