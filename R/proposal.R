#' Randomly and adaptively proposes a new value for the parameter
#'
#' The new value is proposed using a symmetric normal probability density function.
#' Bounded parameters are transformed before running the proposal to cover the full range of R.
#' Similarly, immediately after proposal, the proposed value is transformed back to its
#' original scale to be used in the rest of algorithm. It adaptively searches for sd that result in
#' optimal acceptance rate of 0.44 for 1 dimensional update.
#'
#' @param param      last update of the parameter
#' @param nRun       number of MCMC iterations so far
#' @param accepRate  rate of acceptance so far
#' @param sdLast     the last proposal sd used
#'
#' @return a list of consisting the proposed value of the parameter and the proposal sd used
proposal <- function(param, nRun, accepRate, sdLast) {

  if (nRun <= 5)
    sdNew <- 0.1
  else if (accepRate > 0.44)
    sdNew <- exp(log(sdLast) + min(0.01, sqrt(1/nRun)))
  else
    sdNew <- exp(log(sdLast) - min(0.01, sqrt(1/nRun)))

  if (index %in% cache$ikappa) {
    proposed <- 0.95* rnorm(1, mean = param, sd = sdNew) + 0.05*rnorm(1, mean=param, sd=0.1)

  } else if (index %in% c(cache$ithetaS, cache$ithetaB)) {
    proposed <- exp(0.95* rnorm(1, mean = log(param), sd = sdNew) + 0.05*rnorm(1, mean=log(param), sd=0.1))

  } else if (index %in% c(cache$ialphaS, cache$ialphaB)) {
    proposed <- 1 + (1 / (1 + exp(-(0.95* rnorm(1, mean = log(param-1) - log(2-param), sd = sdNew) + 0.05*rnorm(1, mean=log(param-1) - log(2-param), sd=0.1)))))

  } else if (index %in% c(cache$isigma2S, cache$isigma2B, cache$isigma2E)) {
    proposed <- exp(0.95* rnorm(1, mean = log(param), sd = sdNew) + 0.05*rnorm(1, mean=log(param), sd=0.1))
  }

  return(list(proposed = proposed, sd = sdNew))
}
