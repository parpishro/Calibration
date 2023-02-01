#' setup_prior: based on a given prior type, parameters, branches out to builds
#'              up and return a prior function and its distribution mean (to be
#'              used as initial value) in a list
#'
#' @param pr (String)  type of the prior to be used
#' @param p1 (numeric) first parameter of the prior distribution
#' @param p2 (numeric) second parameter of the prior distribution
#'
#' @return (list(function, numeric)) a list with two member:
#'          - fun : density function of chosen prior distribution
#'          - mean: mean of the chosen distribution (used as initial value)
setup_prior <- function(pr, p1, p2) {
  PRIORS  <- c('Uniform', 'Inverse Gamma', 'Gaussian', 'Beta', 'Gamma',
                'Exponential', 'Logistic', 'Log Beta', 'Jefferys')
  prior   <- match(pr, PRIORS, nomatch = -1)

  if (is.na(prior))
    stop("prior missing!")

  if (prior == -1)
    stop("prior invalid!")

  if (prior == 'Uniform') {
    fun  <- function(x) dunif(x, min = p1, max =  p2, log = TRUE)
    mean <- (p1 + p2) / 2
  }

  if (prior == 'Inverse Gamma') {  # TODO
    fun  <- function(x) 0
    mean <- 0
  }

  if (prior == 'Gaussian') {      # TODO
    fun  <- function(x) 0
    mean <- 0
  }

  if (prior == 'Beta') {          # TODO
    fun  <- function(x) 0
    mean <- 0
  }

  if (prior == 'Gamma') {         # TODO
    fun  <- function(x) 0
    mean <- 0
  }

  if (prior == 'Exponential') {   # TODO
    fun  <- function(x) 0
    mean <- 0
  }

  if (prior == 'Log Beta') {      # TODO
    fun  <- function(x) 0
    mean <- 0
  }

  if (prior == 'Jefferys') {      # TODO
    fun  <- function(x) 0
    mean <- 0
  }

  return(list(fun = fun, mean = mean))
}
