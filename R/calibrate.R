# Calibrates the simulator using both simulator and field data and returns the
# MCMCM distribution of calibration parameters and hyperparameters of the model


# REQUIRE: simulation data (s) must be a matrix with rows representing
#             simulator observations and first column being the univariate
#             response variable, followed by experimental input columns,
#             followed by calibration input columns.
#          field data (field) must be a matrix with rows representing
#             field observations and first column being the univariate response
#             variable, followed by experimental input columns.
#          prior for calibration parameter (theta_pr) is assumed to be "uniform"
#          prior for lambda (lambda_pr), the scale parameter of correlation
#             functions (omega) is assumed to be "logistic"
#          prior for gamma (gamma_pr), the transformed smoothness parameter in
#             correlation functions (alpha), is assumed to be "uniform"
#
#          prior for uncertainty (sigma2_pr) assumed to be "inverse gamma"
# EFFECT: given simulation data, field data, priors for calibration parameters,
#             correlation functions hyperparameters, uncertainty hyperparameters
#             , and MCMC specifications, build a KOH model and run MH MCMC to
#             find the posterior distribution of parameters, hyperparameters,
#             and and their point estimates and estimated variances.
#             The function specifies the data and passes its environment to all
#             child functions for use without needing to pass them as argumnets.
# TODO: 1 - change initail values to vector
#       2 - code few more options for priors

#' calibrate
#'
#' @param sim       (m * (1+p+q) matrix) simulation data
#' @param field     (n * (1+p+) matrix) field data
#' @param Nmcmc     number of MCMC runs
#' @param nBurn     number of MCMC burn ins
#' @param thining   thining rate to de-correlate MCMC results
#' @param theta_pr  prior function for calibration parameters #TODO
#' @param lambda_pr  prior function for scale parameters #TODO
#' @param gamma_pr  prior function for smoothness parameters #TODO
#' @param sigma2_pr prior function for variance parameters #TODO
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
                      Nmcmc = 100, nBurn = 40, thining = 1,
                      theta_pr = "uniform", lambda_pr = "logbeta",
                      gamma_pr = "logistic", sigma2_pr = "inverse gamma") {

  env      <- environment()

  m         <- nrow(sim)               # number of simulation runs
  n         <- nrow(field)             # number of field observations
  p         <- ncol(field) - 1         # number of experimental variables
  q         <- ncol(sim) - p - 1       # number of calibration parameters
  d         <- ncol(sim) - 1           # number of all variables for simulation
  k         <- q + (p + q) + (p + q) +  p + p + 1 + 1 + 1 + 1

  # indices for parameters in phi
  calib     <- 1:q
  scaleS    <- (q+1): (q + (p + q))
  smoothS   <- (q + (p + q) + 1): (q + (p + q) + (p + q))
  scaleB    <- (q + (p + q) + (p + q) + 1): (q + (p + q) + (p + q) + p)
  smoothB   <- (q + (p + q) + (p + q) + p + 1): (k - 4)
  sig2S     <- k - 3
  sig2B     <- k - 2
  sig2E     <- k - 1
  muHat     <- k #

  # number of total parameters =
  #       calibration + sim scale + sim smoothness + bias scale +
  #       bias smoothness + sim variance + bias variance + error variance



  Xs        <- sim[, 1:d]
  ys        <- sim[, d + 1]
  Xb        <- field[, 1:p]
  yf        <- field[, P + 1]
  y         <- (y - mean(ys)) / sd(ys)

  Phi       <- matrix(nrow = Nmcmc, ncol = k)
  Phi[1, ]  <- initialize_phi(env)

  params    <- mcmc(Nmcmc, nBurn, thining, environment = env)

  paramMean <- apply(params, 2, mean)
  paramVar  <- apply(params, 2, var) + apply(sigma_hat, 2, mean)

  results   <- list(mean = paramMean, var = paramVar, distribution = params)
}


