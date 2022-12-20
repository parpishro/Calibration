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
mcmc <- function(ll_fun, prop_fun, Nmcmc, burnIn, phiInit) {
  k <- length(phiInit)
  phi <- matrix(nrow = Nmcmc, ncol = k)
  phi[1,] <- phi_init
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
