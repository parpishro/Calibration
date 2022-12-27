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
#' @param ll_phi log likelihood function of \phi (parameter vector) given data
#' @param prop a proposal function to draw new parameter estimates
#' @param Nmcmc an integer representing total number of MCMC runs
#' @param n_burn an integer representing number of burn-in runs
#' @param phi_init a list of doubles representing initial values for all
#'                 parameters
#'
#' @return a clb object that includes a vector of all estimated parameters
#' @export
#'
#' @examples
mcmc <- function(sim, field,
                 Nmcmc, nBurn, thining,
                 theta_pr, omega_pr, alpha_pr, sigma2_pr,
                 theta0, omega0, alpha0, sigma20) {
  p <- ncol(field) - 1         # number of experimental variables
  q <- ncol(sim) - p           # number of calibration parameters
  # number of total parameters
  k <- q +          # number of calibration parameters
      (p + q) +     # number of scale parameters (for simulation data)
      (p + q) +     # number of smoothness parameters (for simulation data)
       p +          # number of scale parameters (for field data)
       p +          # number of smoothness parameters (for field data)
       1 +          # marginal variance of simulator
       1 +          # marginal variance of bias correction
       1 +          # variance of measurement
       1            # estimate of mean response
  phi      <- matrix(nrow = Nmcmc, ncol = k)
  phi[1,]  <- initialize_phi(k, p, q, theta0, omega0, alpha0, sigma20)
  llPos    <- vector("double", Nmcmc)
  llPos[1] <- log_lik(phi[1,])  # TO DO
  for (i in 2:Nmcmc) {
    for (j in 1:k) {
      param    <- proposal(phi[i-1, j])
      llPos[i] <- update_cov(c(phi[i, 1:j-1],param, phi[i-1, j+1:k]), j)
      u <- unif(1)
      if (llPos[i] - llPos[i - 1] > log(u)) {
        phi[i, j] <- param
      } else {
        phi[i, j] <- phi[i - 1, j]
      }
    }
  }
  return(list(llPos[burnIn:Nmcmc], phi[burnIn:Nmcmc, ]))
}
