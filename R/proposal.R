#' Adaptive Proposal
#'
#'
#' Given a current parameter value, `proposal()` proposes a new value using a symmetric
#' normal probability density function. The mean of the distribution is the current value
#' and its SD is modified adaptively. It modifies SD such that acceptance ratio for the
#' proposed parameter hovers around 0.44 which leads to fastest convergence for
#' one-dimensional updates. Bounded parameters are transformed before running the proposal
#' to expand the possible range to all real numbers. Similarly, immediately after proposal,
#' the proposed value is transformed back to its original scale to be used in the rest of
#' algorithm.
#'
#' @param param      last update of the parameter
#' @param iparam     index of parameter in parameter vector
#' @param iteration  number of MCMC iterations so far
#' @param accepRate  rate of acceptance so far
#' @param sdLast     the last proposal sd used
#'
#' @noRd
#' @return a list consisting of the proposed value of the parameter and the proposal sd used
proposal <- function(param, iparam, iteration, accepRate, sdLast) {
  init <- cache$sdRates[iparam]
  if (iteration <= 5)
    sdNew <- init
  else if (iteration %% 10 == 6 && accepRate > 0.44)
    sdNew <- exp(log(sdLast) + (init/sqrt(iteration)))
  else if (iteration %% 10 == 6 && accepRate < 0.44)
    sdNew <- exp(log(sdLast) - (init/sqrt(iteration)))
  else
    sdNew <- sdLast


  if (iparam %in% c(cache$ikappa, cache$imuB))
    proposed   <- rnorm(1, mean = param, sd = sdNew)

  else if (iparam %in% c(cache$ithetaS, cache$ithetaB, cache$isigma2S, cache$isigma2B, cache$isigma2E))
    proposed   <- exp(rnorm(1, mean = log(param), sd = sdNew))

  else if (iparam %in% c(cache$ialphaS, cache$ialphaB))
    proposed   <- 1 + (1/(1+exp(-rnorm(1, mean = log(param-1) - log(2-param), sd = sdNew))))
  return(list(proposed = proposed, sd = sdNew))
}

