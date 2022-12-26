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
# mcmc <- function(ll_fun, prop_fun, Nmcmc, burnIn, phiInit) {
  p <- ncol(field) - 1         # number of experimental variables
  q <- ncol(sim) - ncol(field) # number of calibration parameters
  k <- 1 + q + (2 * (p + q)) + (2 * p) + 1 + 1 + 1  # number of total parameters
  phi <- matrix(nrow = Nmcmc, ncol = k)
  phi[1,] <- initialize_phi(k, p, q, theta0, omega0, alpha0, sigma20)
  ll <- vector("double", Nmcmc)
  ll[1] <- ll_fun(phi[1,])  # TO DO
  for (i in 2:Nmcmc) {
    for (j in 1:k) {
      phi_ij <- prop_fun(phi[i-1, j])
      ll[i] <- ll_fun(c(phi[i, 1:j-1],phi_ij, phi[i-1, j+1:k]))
      u <- unif(1)
      if (ll[i] - ll[i - 1] > log(u)) {
        phi[i, j] <- phi_ij
      } else {
        phi[i, j] <- phi[i - 1, j]
      }
    }
  }
  return(list(ll[burnIn:Nmcmc], phi[burnIn:Nmcmc, ]))
}
