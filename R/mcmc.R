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
                 thetaPr, lambdaPr, gammaPr, sigma2Pr) {


  # indices for parameters in phi
  itheta     <- get('itheta',   envir = cache)
  ilambdaS   <- get('ilambdaS',  envir = cache)
  igammaS    <- get('igammaS',  envir = cache)
  ilambdaB   <- get('ilambdaB',  envir = cache)
  igammaB    <- get('igammaB',  envir = cache)
  isigma2S   <- get('isigma2S', envir = cache)
  isigma2B   <- get('isigma2B', envir = cache)
  isigma2E   <- get('isigma2E', envir = cache)
  imuHat     <- get('imuHat',   envir = cache)

  # parameters (initialize first row of Phi matrix)
  Phi        <- matrix(nrow = Nmcmc, ncol = cache$k)
  Phi[1, ]   <- init$phi
  logPost    <- double(Nmcmc)
  logPost[1] <- init$logPost

  for (i in 2:Nmcmc) {
    lPost <- logPost[i-1]
    for (j in 1:(cache$k-1)) {
      changed <- proposal(Phi[1:(i-1) ,j])
      if (j == 1)
        params  <- c(changed, Phi[i-1, 2:cache$k])
      else
        params  <- c(Phi[i, 1:j-1], changed, Phi[i-1, (j+1):cache$k])
      chol    <- update_cov(params, j)
      if (is.null(chol)) {
        lPost <- -.Machine$double.xmax
      } else {
        lPost   <- sum(sapply(params[itheta],                thetaPr$fun)  +
                       sapply(params[c(ilambdaS, ilambdaB)], lambdaPr$fun)  +
                       sapply(params[c(igammaS, igammaB)],   gammaPr$fun)  +
                       sapply(params[isigma2S:isigma2E],     sigma2Pr$fun)) -
          ((chol$logDetCov - drop(cache$res %*% chol$InvCov %*% cache$res))/2)
      }

      if (lPost - logPost[i-1] > log(runif(1))) {
        Phi[i, j]  <- changed
        logPost[i] <- lPost
      } else {
        Phi[i, j]  <- Phi[i - 1, j]
        logPost[i] <- logPost[i-1]
      }
    }

  }

  indices <- seq(nBurn, Nmcmc, by = thining)
  Params  <- Phi[indices, ]

  return(list(Params = Params, logPost = logPost))
}
