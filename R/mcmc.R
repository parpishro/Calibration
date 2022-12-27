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
#' @param phiInit  a double vector of initial values for all parameters
#' @param env      an environment of the calibrate function to be set as parent
#'
#' @return a KOH object that includes a matrix of all parameters' distribution
#'            and a vector of log likelihood updates of each usable MCMC runs
mcmc <- function(Nmcmc, nBurn, thining, phiInit, env) {

  env0             <- environment()
  parent.env(env0) <- env

  phi              <- matrix(nrow = Nmcmc, ncol = k)
  phi[1,]          <- phiInit

  logLik           <- vector("double", Nmcmc)
  logLik[1]        <- log_lik(phi[1,], env)


  for (i in 2:Nmcmc) {
    for (j in 1:k) {
      param        <- proposal(phi[i-1, j])
      logLik[i]    <- update_cov(c(phi[i, 1:j-1],param, phi[i-1, j+1:k]), j)

      if (logLik[i] - logLik[i - 1] > log(unif(1))) {
        phi[i, j]  <- param

      } else {
        phi[i, j]  <- phi[i - 1, j]
      }
    }
  }
  indices          <- seq(burnIn:Nmcmc, by = thinning)

  return(list(logLik = logLik[indices], params = phi[indices, ]))
}
