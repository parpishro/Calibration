#' Randomly and adaptively proposes a new value for the parameter
#'
#' The new value is proposed using a symmetric normal probability density function.
#' Bounded parameters are transformed before running the proposal to cover the full range of R.
#' Similarly, immediately after proposal, the proposed value is transformed back to its
#' original scale to be used in the rest of algorithm.
#'
#' @param param sample of parameter so faar
#' @param index index of the parameter
#'
#' @return the proposed value for the parameter
proposal <- function(param, index) {

  if (index %in% cache$ikappa) {
    last     <- param[length(param)]
    sdSoFar  <- if (length(param) > 1) max(0.1, sd(param)) else param
    proposed <- 0.95 * rnorm(1, mean = last, sd = 2.38*sdSoFar)
              + 0.05 * rnorm(1, mean = last, sd = 0.1)

  } else if (index %in% c(cache$ithetaS, cache$ithetaB)) {
    lambda   <- log(param)
    last     <- lambda[length(lambda)]
    sdSoFar  <- if (length(lambda) > 1) max(0.1, sd(lambda)) else lambda
    proposed <- exp(0.95 * rnorm(1, mean = last, sd = 2.38*sdSoFar)
                  + 0.05 * rnorm(1, mean = last, sd = 0.1))

  } else if (index %in% c(cache$ialphaS, cache$ialphaB)) {
    tau      <- log(param-1) - log(2-param)
    last     <- tau[length(tau)]
    sdSoFar  <- if (length(tau) > 1) max(0.1, sd(tau)) else tau
    proposed <- 1 + (1 / (1 + exp(-(0.95 * rnorm(1, mean = last, sd = 2.38*sdSoFar)
                                  + 0.05 * rnorm(1, mean = last, sd = 0.1)))))

  } else if (index %in% c(cache$isigma2S, cache$isigma2B, cache$isigma2B)) {
    logS2    <- log(param)
    last     <- logS2[length(logS2)]
    sdSoFar  <- if (length(logS2) > 1) max(0.1, sd(logS2)) else logS2
    proposed <- exp(0.95 * rnorm(1, mean = last, sd = 2.38*sdSoFar)
                  + 0.05 * rnorm(1, mean = last, sd = 0.1))
  }

  return(proposed)
}
