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
#' @noRd
#' @return a list of consisting the proposed value of the parameter and the proposal sd used
proposal <- function(param, nRun, accepRate, sdLast, index) {



  if (index %in% cache$ikappa) {
    sdNew      <- compute_sd(0.2, nRun, accepRate, sdLast)
    paramTr    <- param
    proposedTr <- rnorm(1, mean = paramTr, sd = sdNew)
    proposed   <- proposedTr
    if(is.na(proposedTr)) stop(paste(nRun, index, sdLast, sdNew, param, paramTr, proposed))

  } else if (index %in% cache$ithetaS) {
    sdNew      <- compute_sd(0.3, nRun, accepRate, sdLast)
    paramTr    <- log(param)
    proposedTr <- rnorm(1, mean = paramTr, sd = sdNew)
    proposed   <- exp(proposedTr)
    if(is.na(proposedTr)) stop(paste(nRun, index, sdLast, sdNew, param, paramTr, proposed))

  } else if (index %in% cache$ithetaB) {
    sdNew      <- compute_sd(0.4, nRun, accepRate, sdLast)
    paramTr    <- log(param)
    proposedTr <- rnorm(1, mean = paramTr, sd = sdNew)
    proposed   <- exp(proposedTr)
    if(is.na(proposedTr)) stop(paste(nRun, index, sdLast, sdNew, param, paramTr, proposed))

  } else if (index %in% cache$ialphaS) {
    sdNew      <- compute_sd(0.3, nRun, accepRate, sdLast)
    paramTr    <- log(param-1) - log(2-param)
    proposedTr <- rnorm(1, mean = paramTr, sd = sdNew)
    proposed   <- 1 + (1 / (1 + exp(-proposedTr)))
    if(is.na(proposedTr)) stop(paste(nRun, index, sdLast, sdNew, param, paramTr, proposed))

  } else if (index %in% cache$ialphaB) {
    sdNew      <- compute_sd(0.4, nRun, accepRate, sdLast)
    paramTr    <- log(param-1) - log(2-param)
    proposedTr <- rnorm(1, mean = paramTr, sd = sdNew)
    proposed   <- 1 + (1 / (1 + exp(-proposedTr)))
    if(is.na(proposedTr)) stop(paste(nRun, index, sdLast, sdNew, param, paramTr, proposed))

  } else if (index == cache$isigma2S) {
    sdNew      <- compute_sd(0.25, nRun, accepRate, sdLast)
    paramTr    <- log(param)
    proposedTr <- rnorm(1, mean = paramTr, sd = sdNew)
    proposed   <- exp(proposedTr)
    if(is.na(proposedTr)) stop(paste(nRun, index, sdLast, sdNew, param, paramTr, proposed))

  } else if (index == cache$isigma2B) {
    sdNew    <- compute_sd(1.5, nRun, accepRate, sdLast)
    paramTr  <- log(param)
    proposedTr <- rnorm(1, mean = paramTr, sd = sdNew)
    proposed   <- exp(proposedTr)
    if(is.na(proposedTr)) stop(paste(nRun, index, sdLast, sdNew, param, paramTr, proposed))

  } else if (index == cache$isigma2E) {
    sdNew    <- compute_sd(0.3, nRun, accepRate, sdLast)
    paramTr  <- log(param)
    proposedTr <- rnorm(1, mean = paramTr, sd = sdNew)
    proposed   <- exp(proposedTr)
    if(is.na(proposedTr)) stop(paste(nRun, index, sdLast, sdNew, param, paramTr, proposed))
  }

  return(list(proposed = proposed, sd = sdNew, proposedTr = proposedTr))
}


compute_sd <- function(init, nRun, accepRate, sdLast) {
  if (nRun <= 5)
    sdNew <- init
  else if (nRun %% 20 == 6 && accepRate > 0.44)
    sdNew <- exp(log(sdLast) + (init/sqrt(nRun)))
  else if (nRun %% 20 == 6 && accepRate < 0.44)
    sdNew <- exp(log(sdLast) - (init/sqrt(nRun)))
  else
    sdNew <- sdLast
  return(sdNew)
}
