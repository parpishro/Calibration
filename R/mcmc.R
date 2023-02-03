# run MH MCMC that computes loglikelihood in each run and updates the model
# parameters, which will lead to a distribution of unkown parameters

# REQUIRE: number of MCMC runs (Nmcmc) must be minimum integer 10000 (default)
#          number of MCMC burns (nBurn) must be an integer minimum 500 (default)
#             and maximum (Nmcmc / 2).
# EFFECT: given simulation data, field data, priors for calibration parameters,
#             correlation functions hyperparameters, uncertainty hyperparameters
#             , and MCMC specifications, run MH MCMC algorithm that compute
#             joint log likelihood in each run to update the parameters. The
#             distribution of each parameter would be the chain of updated
#             parameters. The nBurn results from the start would be omitted to
#             remove the effect of initial start point, and the rest would be
#             thinned by thinning factor to de-correlate the results. In the end
#             each column of the results matrix will represent the posterior
#             distribution of parameters, hyperparameters

#' MH MCMC method for full Bayesian parameter estimation
#'
#' @param Nmcmc    an integer representing total number of MCMC runs
#' @param nBurn    an integer representing number of burn-in runs
#' @param thining  an integer representing number of burn-in runs
#'
#' @return a KOH object that includes a matrix of all parameters' distribution
#'            and a vector of log likelihood updates of each usable MCMC runs
mcmc <- function(Nmcmc, nBurn, thining, init,
                 thetaPr, omegaPr, alphaPr, sigma2Pr) {


  # indices for parameters in phi
  iTheta     <- get('iTheta',   envir = cache)
  iOmegaS    <- get('iOmegaS',  envir = cache)
  iAlphaS    <- get('iAlphaS',  envir = cache)
  iOmegaB    <- get('iOmegaB',  envir = cache)
  iAlphaB    <- get('iAlphaB',  envir = cache)
  iSigma2S   <- get('iSigma2S', envir = cache)
  iSigma2B   <- get('iSigma2B', envir = cache)
  iSigma2E   <- get('iSigma2E', envir = cache)
  iMuHat     <- get('iMuHat',   envir = cache)

  # parameters (initialize first row of Phi matrix)
  Phi        <- matrix(nrow = Nmcmc, ncol = k)
  Phi[1]     <- init$phi
  logPost    <- double(Nmcmc)
  logPost[1] <- init$logPost

  for (i in 2:Nmcmc) {
    lPost <- logPost[i-1]
    for (j in 1:k) {
      changed <- proposal(Phi[1:(i-1) ,j])
      params  <- c(Phi[i, 1:j-1], changed, Phi[i-1, min(k, (j+1)):k])
      chol    <- update_cov(params, j)
      lPost   <- sum(sapply(params[iTheta],              thetaPr$fun)  +
                     sapply(params[c(iOmegaS, iOmegaB)], omegaPr$fun)  +
                     sapply(params[c(iAlphaS, iAlphaB)], alphaPr$fun)  +
                     sapply(params[iSigma2S:iSigma2E],   sigma2Pr$fun)) -
                ((chol$logDetCov - (cache$res %*% chol$InvCov %*% cache$res))/2)

      if ((lPost - logPost[i-1]) > log(runif(1)))
        Phi[i, j]  <- changed
      else
        Phi[i, j]  <- Phi[i - 1, j]
    }
    logPost[i] <- lPost
  }

  indices <- seq(burnIn:Nmcmc, by = thinning)
  Params  <- Phi[indices, ]

  return(list(Params = Params, logPost = logPost))
}
