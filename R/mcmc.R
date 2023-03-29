# Runs an adaptive Metropolis-within-Gibbs MCMC that computes posterior log
# likelihood in each run and updates the model parameters. New parameters are
#' proposed according to a symmetric and adaptive proposal function. After
#' removing burn-in runs and thinning the samples, the output approximates the
#' posterior distribution of parameters.

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
#' @param Nmcmc     an integer representing total number of MCMC runs
#' @param nBurn     an integer representing number of burn-in runs
#' @param thining   an integer representing number of burn-in runs
#' @param init      a list containing initial values for parameters (created by setup_prior), the first update of log posterior
#' @param kappaPr   a function that computes log prior for calibration parameters
#' @param thetaPr   a function that computes log prior for cscale parameters
#' @param alphaPr   a function that computes log prior for smoothness parameters
#' @param sigma2Pr  a function that computes log prior for calibration parameters
#' @noRd
#' @return a KOH object that includes a matrix of all parameters' distribution
#'            and a vector of log likelihood updates of each usable MCMC runs
mcmc <- function(Nmcmc, nBurn, thining, init,
                 kappaPr, thetaPr, alphaPr, sigma2Pr) {

  # indices for parameters in phi
  ikappa     <- get('ikappa',   envir = cache)
  ithetaS    <- get('ithetaS',  envir = cache)
  ialphaS    <- get('ialphaS',  envir = cache)
  ithetaB    <- get('ithetaB',  envir = cache)
  ialphaB    <- get('ialphaB',  envir = cache)
  isigma2S   <- get('isigma2S', envir = cache)
  isigma2B   <- get('isigma2B', envir = cache)
  isigma2E   <- get('isigma2E', envir = cache)
  #imuHat     <- get('imuHat',   envir = cache)

  # parameters (initialize first row of Phi matrix)
  Phi        <- matrix(nrow = Nmcmc, ncol = cache$k)
  Phi[1, ]   <- init$phi
  cat("initial values: ", round(init$phi, 3), "\n")
  logPost    <- double(Nmcmc)
  logPost[1] <- init$logPost

  indices    <- seq(nBurn, Nmcmc, by = thining)

  accepLast  <- double(cache$k)
  quantAccep <- double(cache$k)
  accepSoFar <- double(cache$k)
  sdProp     <- double(cache$k)
  for (i in 2:Nmcmc) {
    logPost[i] <- logPost[i-1]
    lPost      <- logPost[i-1]
    for (j in 1:(cache$k)) {
      res <- proposal(Phi[i-1 ,j], i, quantAccep[j], sdProp[j], j)
      changed   <- res$proposed
      sdProp[j] <- res$sd
      if (j == 1)
        params <- c(changed, Phi[i-1, 2:cache$k])
      else if (j == cache$k)
        params <- c(Phi[i, 1:(cache$k-1)], changed)
      else
        params <- c(Phi[i, 1:j-1], changed, Phi[i-1, (j+1):cache$k])

      chol <- update_cov(params, iChanged = j)

      if (is.null(chol))
        lPost <- -.Machine$double.xmax
      else
        lPost <- sum(sapply(params[ikappa],              kappaPr)  +
                     sapply(params[c(ithetaS, ithetaB)], thetaPr)  +
                     sapply(params[c(ialphaS, ialphaB)], alphaPr)  +
                     sapply(params[isigma2S:isigma2E],   sigma2Pr)) -
                 ((chol$logDetCov - drop(t(cache$y) %*% chol$InvCov %*% cache$y))/2)

      if (is.finite(lPost) & lPost - logPost[i] >  log(runif(1))) {
        Phi[i, j]     <- changed
        logPost[i]    <- lPost
        accepSoFar[j]  <- accepSoFar[j] + 1
      } else {
        Phi[i, j]     <- Phi[i - 1, j]
      }
    }
    #Phi[i, imuHat] <- update_mu()

    if (i %% floor(Nmcmc/1000) == 0 && i/floor(Nmcmc/1000) <= 1000) {   #
      cat("------------------------------------------------------------------------", "\n")
      cat("finished ",  (i/floor(Nmcmc/1000))/10, "% of MCMC runs...", "\n")
      quantAccep <- (accepSoFar - accepLast)*1000/Nmcmc
      cat("acceptance ratio in the last percent: ", round(quantAccep, 2), "\n")
      cat("parameter sample:                     ", round(Phi[i, ], 4), "\n")
      cat("proposal sd:                          ", round(sdProp, 4), "\n")
      cat("total acceptance ratio:               ", round(accepSoFar/i, 2), "\n")
      accepLast <- accepSoFar
    }
  }
  cat("Completed MCMC runs.", "\n")
  cat("total acceptance ratio ", round(accepSoFar/Nmcmc, 3), "\n")
  cat("parameter sample: ", round(Phi[Nmcmc, ], 3), "\n")
  assign('Params',     Phi[indices, ],   envir = cache)
  assign('Phi',        Phi,              envir = cache)
  assign('logPost',    logPost[indices], envir = cache)
  assign('acceptance', accepSoFar,       envir = cache)
}
