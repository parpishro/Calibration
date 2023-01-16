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
#' @param env      an environment of the calibrate function to be set as parent
#'
#' @return a KOH object that includes a matrix of all parameters' distribution
#'            and a vector of log likelihood updates of each usable MCMC runs
mcmc <- function(Nmcmc, nBurn, thining, env) {

  env0             <- environment()
  parent.env(env0) <- env

  Xf    <- cbind(Xb, replicate(n, Phi[1, calib], n))
  CorFF <- corelation(Xf, Xf, scale = Phi[1, scaleS], smooth = Phi[1, smoothS])
  CorFS <- corelation(Xf, Xs, scale = Phi[1, scaleS], smooth = Phi[1, smoothS])
  CorSF <- t(CorFS)
  CorSS <- corelation(Xs, Xs, scale = Phi[1, scaleS], smooth = Phi[1, smoothS])
  CorB  <- corelation(Xb, Xb, scale = Phi[1, scaleB], smooth = Phi[1, smoothB])
  muHat <- mu_hat(CorSS, ys)
  res   <- y - muHat



  sig2S <- Phi[1, sig2S]
  sig2B <- Phi[1, sig2B]
  sig2E <- Phi[1, sig2E]


  I         <- diag(n)
  AugCov    <- rbind(cbind(((sig2S * CorFF) + (sig2B * CorB) + (sig2E * I)), (sig2S * CorFS)),
                     cbind((sig2S * CorSF), (sig2S * CorSS)))
  logPost   <- double(((Nmcmc - 1) * k) + 1)
  logPost[1]<- log_prior(Phi[1, calib],
                         Phi[1, c(scaleS, scaleB)],
                         Phi[1, c(smoothS, smoothB)],
                         Phi[1, c(sig2S, sig2B, sig2E)])
                - log_lik(chol_cov(env0), res)


  for (i in 2:Nmcmc) {
    for (j in 1:k) {
      print(j)
      changed      <- proposal(Phi[1:(i-1) ,j])
      params       <- c(Phi[i, 1:j-1], changed, Phi[i-1, min(k, (j+1)):k])
      chol         <- update_cov(params, j, env0)
      ind          <- ((i-2) * k) + j + 1
      logPost[ind] <- log_prior(params[calib],
                                params[c(scaleS, scaleB)],
                                params[c(smoothS, smoothB)],
                                params[c(sig2S, sig2B, sig2E)])
                      - log_lik(chol, res)

      if ((logPost[i] - logPost[i - 1]) > log(runif(1))) {
        Phi[i, j]  <- changed

      } else {
        Phi[i, j]  <- Phi[i - 1, j]
      }
    }
  }
  indices <- seq(burnIn:Nmcmc, by = thinning)

  return(list(logPost = logPost[indices], params = Phi[indices, ]))
}
