#' setup_prior: based on a given prior type, parameters, branches out to builds
#'              up and return a log prior function and its distribution mean (to
#'              be used as initial value) in a list
#'
#' @param pr (String)  type of the prior to be used
#' @param p1 (numeric) first parameter of the prior distribution
#' @param p2 (numeric) second parameter of the prior distribution
#'
#' @return (list(function, numeric)) a list with two member:
#'          - fun : density function of chosen prior distribution
#'          - mean: mean of the chosen distribution (used as initial value)
setup_prior <- function(pr, p1, p2) {
  PRIORS  <- c('uniform', 'inversegamma', 'gaussian', 'beta', 'gamma',
                'exponential', 'logistic', 'logbeta', 'jefferys', 'chen')
  prior   <- match(pr, PRIORS, nomatch = -1)
  if (is.na(prior)) {
    stop("prior missing!")
  } else if (prior == -1) {
    stop("prior invalid!")
  } else if (prior == 1) {  # uniform
    # !!! ALL FUNCTIONS RETURN THE LOG VALUE OF THE GIVEN DISTRIBUTION
    fun  <- function(x) dunif(x, min = p1, max =  p2, log = TRUE)
    m    <- (p1 + p2) / 2
  } else if (prior == 2) { # inversegamma
    # REQUIRE: p1 > 0, x > 0
    fun  <- function(x) p1*log(p2) - lgamma(p1) - (p1+1)*log(abs(x)) - p2/(abs(x))
    m    <- p2 / (p1 - 1)
  } else if (prior == 3) { # gaussian
    fun  <- function(x) -log(sqrt(2*p2*pi)) - (((x - p1)^2)/ (2*p2))
    m    <- p1
  } else if (prior == 4) { # beta
    fun  <- function(x) ((p1-1) * log(x)) + ((p2-1) * (1-x)) - lbeta(p1, p2)
    m    <- p1 / (p1+p2)
  } else if (prior == 5) { # gamma
    fun  <- function(x) log(p2^p1) + ((p1-1)*log(x)) - (p2*x) - lgamma(p1)
    m    <- p1 / p2
  } else if (prior == 6) { # exponential
    fun  <- function(x) log(p1) - (p1*x)
    m    <- 1/p2
  } else if (prior == 7) { # logistic    # TODO
    fun  <- function(x) 0
    m    <- 0
  } else if (prior == 8) { # logbeta     # TODO
    fun  <- function(x) 0
    m    <- 0
  } else if (prior == 9) { # jefferys   # TODO
    fun  <- function(x) 0
    m    <- 0
  } else if (prior == 10) { # chen
    # from "Flexible Correlation ..." at 10.1137/15M1008774, P606, Eq. 4.1
    # a log is taken from the equation and simplified
    fun  <- function(x) log(0.5) + 1.5*x + 0.5*log(1+exp(-x)) - 2*log(1+exp(x))
    m    <- 1 # computed numerically
  }

  return(list(fun = fun, mean = m))
}
