# Calibrates the simulator using both simulator and field data and returns the
# distribution of calibration parameters and hyperparameters of the model


# REQUIRE: simulation data (sim) must be a matrix with rows representing
#             simulator observations and first column being the univariate
#             response variable, followed by experimental input columns,
#             followed by calibration input columns.
#          field data (field) must be a matrix with rows representing
#             field observations and first column being the univariate response
#             variable, followed by experimental input columns.
#          prior type for calibration parameters (theta_pr) must be either
#             "weak" (default) or "strong".
#          prior for scale of correlation functions (omega_pr) must be "logbeta"
#          prior for smoothness correlation function (alpha_pr) must be
#             "logistic"
#          prior for uncertainty (sigma2_pr) must be "inverse gamma"
# EFFECT: given simulation data, field data, priors for calibration parameters,
#             correlation functions hyperparameters, uncertainty hyperparameters
#             , and MCMC specifications, build a KOH model and run MH MCMC to
#             find the posterior distribution of parameters, hyperparameters,
#             and and their point estimates and estimated variances.
#             The function specifies the data and passes its environment to all
#             child functions for use without needing to pass them as argumnets.

#' calibrate
#'
#' @param sim       (m * (1+p+q) matrix) simulation data
#' @param field     (n * (1+p+) matrix) field data
#' @param Nmcmc     number of MCMC runs
#' @param nBurn     number of MCMC burn ins
#' @param thining   thining rate to de-correlate MCMC results
#' @param theta_pr  prior function for calibration parameters #TODO
#' @param omega_pr  prior function for scale parameters #TODO
#' @param alpha_pr  prior function for smoothness parameters #TODO
#' @param sigma2_pr prior function for variance parameters #TODO
#' @param theta0    initial value of calibration parameters (to be given MCMC)
#' @param omega0    initial value of scale parameters (to be given MCMC)
#' @param alpha0    initial value of smoothness parameters (to be given MCMC)
#' @param sigma20   initial value of smoothness parameters (to be given MCMC)
#'
#' @return a list containing posterior:
#'    - (((Nmcmc - nBurn) / thinning) * k) parameters distribution:
#'        - 1st q columns: distribution of calibration parameters
#'        - next (p + q) columns: distribution of sim scale parameters
#'        - next (p + q) columns: distribution of sim smoothness parameters
#'        - next p columns: distribution of bias scale parameters
#'        - next p columns: distribution of bias smoothness parameters
#'        - next column (index: k - 3): distribution of sim variance parameter
#'        - next column (index: k - 1): distribution of bias variance parameter
#'        - next column (index: k - 2): distribution of error variance parameter
#'        - next column (index: k): distribution of mean response parameter
#'    - a vector (length: k) of parameter estimated variances
#'    - a vector (length: k) of parameter point estimates
#'
#' @export
#'
#' @examples
calibrate <- function(sim, field,
                      Nmcmc = 10000, nBurn = 500, thining = 100,
                      theta_pr = "weak", omega_pr = "logbeta",
                      alpha_pr = "logistic", sigma2_pr = "inverse gamma",
                      theta0, omega0, alpha0, sigma20) {

  env      <- environment()

  m         <- nrow(sim)               # number of simulation runs
  n         <- nrow(field)             # number of field observations
  p         <- ncol(field) - 1         # number of experimental variables
  q         <- ncol(sim) - p - 1       # number of calibration parameters
  d         <- ncol(sim) - 1           # number of all variables for simulation

  # number of total parameters =
  #       calibration + sim scale + sim smoothness + bias scale +
  #       bias smoothness + sim variance + bias variance + error variance
  k         <- q + (p + q) + (p + q) +  p + p + 1 + 1 + 1 + 1
  phiInit   <- initialize_phi(k, p, q, theta0, omega0, alpha0, sigma20)

  Xs        <- sim[, 1:d]
  ys        <- sim[, d + 1]
  Xb        <- field[, 1:p]
  yf        <- field[, P + 1]
  y         <- c(ys, yf)

  theta_pr  <- theta_pr
  omega_pr  <- omega_pr
  alpha_pr  <- alpha_pr
  sigma2_pr <- sigma2_pr

  params    <- mcmc(Nmcmc, nBurn, thining, phiInit, environment = env)
  mu-hat    <- clb$params[, k]
  sigma_hat <- clb$params[, (k-1)]
  paramVar  <- apply(mu-hat, 2, var) + apply(sigma_hat, 2, mean)
  paramMean <- apply(mu-hat, 2, mean)
  results   <- list(mean = paramMean, var = paramVar)
}


