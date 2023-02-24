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
#             functions (lambda) is assumed to be "logistic"
#          prior for gamma (gamma_pr), the transformed smoothness parameter in
#             correlation functions (gamma), is assumed to be "uniform"
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
#   MCMCM distribution of calibration parameters and hyperparameters of the model
#'
#' @param sim       (m * (1+p+q) matrix) simulation data
#' @param field     (n * (1+p+) matrix) field data
#' @param Nmcmc     number of MCMC runs
#' @param nBurn     number of MCMC burn ins
#' @param thining   thining rate to de-correlate MCMC results
#' @param theta     prior type for calibration parameters
#' @param t1        first parameter of the chosen theta prior distribution
#' @param t2        second parameter of the chosen theta prior distribution
#' @param lambda    prior type for scale parameters
#' @param l1        first parameter of the chosen rho prior distribution
#' @param l2        second parameter of the chosen rho prior distribution
#' @param gamma     prior type for smoothness parameters
#' @param g1        first parameter of the chosen nu prior distribution
#' @param g2        second parameter of the chosen nu prior distribution
#' @param sigma2    prior type for variance parameters
#' @param s1        first parameter of the chosen sigma2 prior distribution
#' @param s2        second parameter of the chosen sigma2 prior distribution
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
                      Nmcmc  = 10000, nBurn = 500, thining = 100,
                      theta  = "uniform",      t1 = 0,   t2 = 20,
                      lambda = "chen",         l1 = NA,  l2 = NA,
                      gamma  = "uniform",      g1 = -20, g2 = 20,
                      sigma2 = "inversegamma", s1 = 2,   s2 = 1) {

  thetaPr    <- setup_prior(theta,  t1, t2)
  # parameters l1, l2 are already factored in 'chen' prior and will not be used
  lambdaPr   <- setup_prior(lambda, l1, l2)
  gammaPr    <- setup_prior(gamma,  g1, g2)
  sigma2Pr   <- setup_prior(sigma2, s1, s2)


  init       <- setup_cache(sim, field, thetaPr, lambdaPr, gammaPr, sigma2Pr)
  Result     <- mcmc(Nmcmc, nBurn, thining, init,
                     thetaPr, lambdaPr, gammaPr, sigma2Pr)

  #paramMean <- apply(Result$Params, 2, mean)
  #paramVar  <- apply(Result$Params, 2, var) + apply(sigma_hat, 2, mean)

  #Resul$distribution[, 2] = 1 - 1/(1+exp(Resul$distribution[, 2]))

  return(list(distributions = Result$Params, logPost = Result$logPost))
}


