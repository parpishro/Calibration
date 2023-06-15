#' Randomly and adaptively proposes a new value for the parameter
#'
#' The new value is proposed using a symmetric normal probability density function.
#' Bounded parameters are transformed before running the proposal to cover the full range of R.
#' Similarly, immediately after proposal, the proposed value is transformed back to its
#' original scale to be used in the rest of algorithm. It adaptively searches for sd that result in
#' optimal acceptance rate of 0.44 for 1 dimensional update.
#'
#' @param param      last update of the parameter
#' @param iteration       number of MCMC iterations so far
#' @param accepRate  rate of acceptance so far
#' @param sdLast     the last proposal sd used
#' @noRd
#' @return a list of consisting the proposed value of the parameter and the proposal sd used
proposal <- function(param, iparam, iteration, accepRate, sdLast) {

  if (iparam %in% cache$ikappa) {
    sdNew      <- compute_sd(0.2, iteration, accepRate, sdLast)
    proposed   <- rnorm(1, mean = param, sd = sdNew)

  } else if (iparam %in% c(cache$ithetaS, cache$ithetaB)) {
    sdNew      <- compute_sd(0.3, iteration, accepRate, sdLast)
    proposed   <- exp(rnorm(1, mean = log(param), sd = sdNew))

  } else if (iparam %in% c(cache$ialphaS, cache$ialphaB)) {
    sdNew      <- compute_sd(0.3, iteration, accepRate, sdLast)
    proposed   <- 1 + (1/(1+exp(-rnorm(1, mean = log(param-1) - log(2-param), sd = sdNew))))

  } else if (iparam == cache$isigma2S) {
    sdNew      <- compute_sd(0.1, iteration, accepRate, sdLast)
    proposed   <- exp(rnorm(1, mean = log(param), sd = sdNew))

  } else if (iparam == cache$isigma2B) {
    sdNew      <- compute_sd(1, iteration, accepRate, sdLast)
    proposed   <- exp(rnorm(1, mean = log(param), sd = sdNew))

  } else if (iparam == cache$isigma2E) {
    sdNew      <- compute_sd(0.1, iteration, accepRate, sdLast)
    proposed   <- exp(rnorm(1, mean = log(param), sd = sdNew))

  } else if (iparam == cache$imu) {
  sdNew      <- compute_sd(1, iteration, accepRate, sdLast)
  proposed   <- rnorm(1, mean = param, sd = sdNew)
  }

  return(list(proposed = proposed, sd = sdNew))
}


compute_sd <- function(init, iteration, accepRate, sdLast) {
  if (iteration <= 5)
    sdNew <- init
  else if (iteration %% 10 == 6 && accepRate > 0.44)
    sdNew <- exp(log(sdLast) + (init/sqrt(iteration)))
  else if (iteration %% 10 == 6 && accepRate < 0.44)
    sdNew <- exp(log(sdLast) - (init/sqrt(iteration)))
  else
    sdNew <- sdLast
  return(sdNew)
}
