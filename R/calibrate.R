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
#' @description  Calibrates the simulator using both simulator and field data and returns the
#  MCMCM distribution of calibration parameters and hyperparameters of the model
#'
#' @param sim       (m * (1+p+q) matrix) simulation data
#' @param field     (n * (1+p+) matrix) field data
#' @param Nmcmc     number of MCMC runs
#' @param nBurn     number of MCMC burn ins
#' @param thining   thining rate to de-correlate MCMC results
#' @param theta_pr  prior function for calibration parameters
#' @param lambda_pr  prior function for scale parameters
#' @param gamma_pr  prior function for smoothness parameters
#' @param sigma2_pr prior function for variance parameters
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
                      Nmcmc  = 100, nBurn = 40, thining = 2,
                      theta  = "uniform", t1, t2,
                      omega  = "logbeta", o1, o2,
                      alpha  = "logistic", a1, a2,
                      sigma2 = "inverse gamma", s1, s2) {

  thetaPr    <- setup_prior(theta,  t1, t2)
  omegaPr    <- setup_prior(omega,  o1, o2)
  alphaPr    <- setup_prior(alpha,  a1, a2)
  sigma2Pr   <- setup_prior(sigma2, s1, s2)

  init       <- setup_cache(sim, field, thetaPr, omegaPr, alphaPr, sigma2Pr)
  params     <- mcmc(Nmcmc, nBurn, thining, init,
                     thetaPr, omegaPr, alphaPr, sigma2Pr)

  paramMean <- apply(params, 2, mean)
  paramVar  <- apply(params, 2, var) + apply(sigma_hat, 2, mean)

  return(list(mean = paramMean, var = paramVar, distribution = params))
}


