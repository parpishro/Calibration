#' Building Prior Function
#'
#' `build_prior` uses user-given prior types and parameters to select and parametrize
#' the prior function. This function
#' setup prior functions to return log of prior distribution.
#'
#' @param dist character string representing prior distribution family
#' @param p1    double representing the first parameter of the prior distribution
#' @param p2    double representing the second parameter of the prior distribution
#'
#' @return      function that given its input, x, computes the log of the chosen prior
#'              probability density function
#' @import      stats
#' @example     man/examples/ex_build_prior.R
#'
#' @export
build_prior <- function(dist, p1, p2) {
  force(p1)
  force(p2)

  if (dist == "uniform") {
    mean <- (p1 + p2)/2
    sd   <- (p2 - p1)/sqrt(12)
    fun  <- function(x) dunif(x, min = p1, max = p2, log = T)

  } else if (dist == "normal") {
    mean <- p1
    sd   <- p2
    fun  <- function(x) dnorm(x, mean = p1, sd = p2, log = T)

  } else if (dist == "normalTr") {
    mean <- p1
    sd   <- p2
    fun  <- function(x) dunif(x, min  = p1 - 2*p2, max = p1 + 2*p2, log = T) +
                        dnorm(x, mean = p1,        sd  = p2,        log = T)

  } else if (dist == "lognormal") {
    mean <- exp(p1 + (p2^2)/2)
    sd   <- sqrt(exp(p2^2) - 1) * exp(p1 + (p2^2)/2)
    fun  <- function(x) dlnorm(x, meanlog  = p1, sdlog  = p2, log = T)

  } else if (dist == "gamma") {
    mean <- p1 / p2
    sd   <- sqrt(p1) / p2
    fun  <- function(x) dgamma(x, shape = p1, scale = p2, log = T)

  } else if (dist == "inversegamma") {
    mean <- if (p1 > 1) p2/(p1 - 1) else 0.1
    sd   <- if (p1 > 2) p2/((p1 - 1)*sqrt(p1 - 2)) else mean/10
    fun  <- function(x) -(p1 + 1)*log(x) - p2/x

  } else if (dist == "beta") {
    mean <- p1/(p1 + p2)
    sd   <- sqrt(p1*p2)/((p1 + p2)*sqrt(p1 + p2 + 1))
    fun  <- function(x) dbeta(x, shape1 = p1, shape2 = p2, log = T)

  } else if (dist == "betashift") {
    mean <- 1 + (p1/(p1 + p2))
    sd   <- sqrt(p1*p2)/((p1 + p2)*sqrt(p1 + p2 + 1))
    fun  <- function(x) dbeta(x - 1, shape1 = p1, shape2 = p2, log = T)

  } else if (dist == "logbeta") {
    mean <- log(p1 + p2) - log(p1)
    sd   <- log(p1 + p2) + 0.5*log(p1 + p2 + 1) - 0.5*(log(p1*p2))
    fun  <- function(x) -p1*x + (p2 - 1)*log(1 - exp(-x))

  } else if (dist == "logistic") {
    mean <- p1
    sd   <- (p2*pi)/sqrt(3)
    fun  <- function(x) dlogis(x, location = p1, scale  = p2, log = T)

  } else if (dist == "exponential") {
    mean <- 1/p1
    sd   <- 1/p1
    fun  <- function(x) dexp(x, rate = p1, log = T)

  } else if (dist == "fixed") {
    mean <- p1
    sd   <- 0
    fun  <- function(x) 0

  } else stop("at least one prior missing or invalid!")


  return(list(mean = mean, sd = sd, fun = fun))
}


