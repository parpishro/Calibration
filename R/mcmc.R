# Runs an adaptive Metropolis-within-Gibbs MCMC that computes posterior log
# likelihood in each run and updates the model parameters. New parameters are
#' proposed according to a symmetric and adaptive proposal function. After
#' removing burn-in runs and thinning the samples, the output approximates the
#' posterior distribution of parameters.
#' # REQUIRE: number of MCMC runs (Nmcmc) must be minimum integer 10000 (default)
#'          number of MCMC burns (nBurn) must be an integer minimum 500 (default)
#'             and maximum (Nmcmc / 2).
#' EFFECT: given simulation data, field data, priors for calibration parameters,
#'             correlation functions hyperparameters, uncertainty hyperparameters
#'             , and MCMC specifications, run MH MCMC algorithm that compute
#'             joint log likelihood in each run to update the parameters. The
#'             distribution of each parameter would be the chain of updated
#'             parameters. The nBurn results from the start would be omitted to
#'             remove the effect of initial start point, and the rest would be
#'             thinned by thinning factor to de-correlate the results. In the end
#'             each column of the results matrix will represent the posterior
#'             distribution of parameters, hyperparameters
#' MH MCMC method for full Bayesian parameter estimation
#'
#' @param Nmcmc     an integer representing total number of MCMC runs
#' @param nBurn     an integer representing number of burn-in runs
#' @param init
#' @param thining   an integer representing number of burn-in runs
#'
#' @noRd
mcmc <- function(Nmcmc, nBurn, thining, init) {
  Phi        <- init$Phi
  logPost    <- init$logPost
  InvCovRes  <- init$InvCovRes
  indices    <- seq(nBurn+1, Nmcmc, by = thining)
  l          <- cache$l

  accepLast  <- double(l)
  accepted   <- double(l)
  accepRate  <- double(l)
  sdProp     <- double(l)
  ones       <- rep(1, cache$m+cache$n)
  res        <- cache$y-Phi[1, cache$imu]


  for (i in 2:Nmcmc) {
    logPost[i]     <- logPost[i-1]
    InvCovRes[i,]  <- InvCovRes[i-1, ]
    for (j in 1:l) {
      if (j %in% cache$ifixed) next
      out         <- proposal(Phi[i-1, j], j, i, accepRate[j], sdProp[j])
      changed     <- out$proposed
      sdProp[j]   <- out$sd

      if (j == 1) {
        params <- c(changed, Phi[i-1, 2:l])
        chol <- update_cov(params, iChanged = j)
      } else if (j == l) {
        params <- c(Phi[i, 1:(l-1)], changed)
        res       <- cache$y-params[cache$imu]
      } else {
        params <- c(Phi[i, 1:(j-1)], changed, Phi[(i-1), (j+1):l])
        chol <- update_cov(params, iChanged = j)
      }


      if (is.null(chol))
        lPost <- -.Machine$double.xmax
      else {
        invCovRes <- chol$InvCov %*% res
        lPost     <- sum(sapply(params[cache$ikappa],  cache$kappa_pr),
                         sapply(params[cache$ithetaS], cache$thetaS_pr),
                         sapply(params[cache$ialphaS], cache$alphaS_pr),
                         sapply(params[cache$ithetaB], cache$thetaB_pr),
                         sapply(params[cache$ialphaB], cache$alphaB_pr),
                         cache$sigma2S_pr(params[cache$isigma2S]),
                         cache$sigma2B_pr(params[cache$isigma2B]),
                         cache$sigma2E_pr(params[cache$isigma2E]),
                         cache$mu_pr(params[cache$imu])) -
                   (0.5*chol$logDetCov) -
                   (0.5*drop(t(res)%*%invCovRes))
      }


      if (is.finite(lPost) && (lPost - logPost[i] >  log(runif(1)))) {
        Phi[i, j]   <- changed
        logPost[i]  <- lPost
        accepted[j] <- accepted[j] + 1
        InvCovRes[i, ] <- invCovRes
      } else {
        Phi[i, j]      <- Phi[i-1, j]
      }
    }

    if (i == 5) {
      accepRate  <- (accepted - accepLast)/4
      accepLast  <- accepted
    }
    if (i > 5 && i %% 10 == 5) {
      accepRate  <- (accepted - accepLast)/10
      accepLast  <- accepted
    }


    if (i %% floor(Nmcmc/100) == 0 && i/floor(Nmcmc/100) <= 100) {   #
      cat("------------------------------------------------------------------------", "\n")
      cat("finished ",  (i/floor(Nmcmc/100)), "% of MCMC runs...", "\n")
      cat("acceptance ratio in the last batch:   ", round(accepRate, 2), "\n")
      cat("parameter sample:                     ", round(Phi[i, ], 4), "\n")
      cat("proposal sd:                          ", round(sdProp, 4), "\n")
      cat("total acceptance ratio:               ", round(accepted/i, 2), "\n")
    }
  }
  cat("------------------------------------------------------------------------", "\n")
  cat("Completed MCMC runs.", "\n")
  cat("total acceptance ratio ", round(accepted/Nmcmc, 3), "\n")
  cache$Params     <- Phi[indices, ]
  cache$logPost    <- logPost[indices]
  cache$acceptance <- accepted/Nmcmc
  cache$InvCovRes  <- InvCovRes[indices, ]

}
