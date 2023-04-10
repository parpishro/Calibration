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
#' @param kappaPr   a function that computes log prior for calibration parameters
#' @param thetaPr   a function that computes log prior for cscale parameters
#' @param alphaPr   a function that computes log prior for smoothness parameters
#' @param sigma2Pr  a function that computes log prior for calibration parameters
#' @noRd
#' @return a KOH object that includes a matrix of all parameters' distribution
#'            and a vector of log likelihood updates of each usable MCMC runs
mcmc <- function(Nmcmc, nBurn, thining,
                 kappaPr, thetaPr, alphaPr, sigma2Pr) {


  indices    <- seq(nBurn, Nmcmc, by = thining)

  accepLast  <- double(cache$k)
  quantAccep <- double(cache$k)
  cache$accepSoFar <- double(cache$k)
  paramTr    <- double(cache$k)
  sdProp     <- double(cache$k)
  for (i in 2:Nmcmc) {
    cache$logPost[i] <- cache$logPost[i-1]
    for (j in 1:(cache$k)) {
      res         <- proposal(cache$Phi[i-1 ,j], i, quantAccep[j], sdProp[j], j)
      changed     <- res$proposed
      paramTr[j]  <- res$proposedTr

      sdProp[j]   <- res$sd

      if (j == 1)
        params <- c(changed, cache$Phi[i-1, 2:cache$k])
      else if (j == cache$k)
        params <- c(cache$Phi[i, 1:(cache$k-1)], changed)
      else
        params <- c(cache$Phi[i, 1:(j-1)], changed, cache$Phi[(i-1), (j+1):cache$k])

      chol <- update_cov(params, iChanged = j)

      if (is.null(chol))
        lPost <- -.Machine$double.xmax
      else
        lPost <- sum(sapply(params[cache$ikappa],              kappaPr)  +
                     sapply(params[c(cache$ithetaS, cache$ithetaB)], thetaPr)  +
                     sapply(params[c(cache$ialphaS, cache$ialphaB)], alphaPr)  +
                     sapply(params[cache$isigma2S:cache$isigma2E],   sigma2Pr)) -
                 ((chol$logDetCov + drop(t(cache$y) %*% chol$InvCov %*% cache$y))/2)

      if (is.finite(lPost) && (lPost - cache$logPost[i] >  log(runif(1)))) {
        cache$Phi[i, j]      <- changed
        cache$logPost[i]     <- lPost
        cache$accepSoFar[j]  <- cache$accepSoFar[j] + 1
      } else {
        cache$Phi[i, j]      <- cache$Phi[i - 1, j]
      }
    }

    if (i == 5) {
      quantAccep <- (cache$accepSoFar - accepLast)/4
      accepLast  <- cache$accepSoFar
    }
    if (i > 5 && i %% 20 == 5) {
      quantAccep <- (cache$accepSoFar - accepLast)/20
      accepLast  <- cache$accepSoFar
    }


    if (i %% floor(Nmcmc/100) == 0 && i/floor(Nmcmc/100) <= 100) {   #
      cat("------------------------------------------------------------------------", "\n")
      cat("finished ",  (i/floor(Nmcmc/100)), "% of MCMC runs...", "\n")
      cat("acceptance ratio in the last batch:   ", round(quantAccep, 2), "\n")
      cat("parameter sample:                     ", round(cache$Phi[i, ], 4), "\n")
      cat("proposal sd:                          ", round(sdProp, 4), "\n")
      cat("total acceptance ratio:               ", round(cache$accepSoFar/i, 2), "\n")
    }
  }
  cat("------------------------------------------------------------------------", "\n")
  cat("Completed MCMC runs.", "\n")
  cat("total acceptance ratio ", round(cache$accepSoFar/Nmcmc, 3), "\n")
  cache$Params <- cache$Phi[indices, ]
}
