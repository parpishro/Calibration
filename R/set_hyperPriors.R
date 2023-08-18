#' Set Calibration Hyperparameter Priors
#'
#' `set_hyperPriors()` can be used to specify one or more hyperparameter priors. There are
#' 8 classes of hyperparameters. Each class can be specified using a named list with four
#' fields:
#'  - `dist`: character string (vector of character strings) to specify the family of the
#'            distributions. All common probability distribution are implemented.
#'  - `p1`:   first distribution parameter
#'  - `p2`:   second distribution parameter
#' All hyperparameters have default values to work with variety of models. Therefore, it
#' is recommended to use this function without any argument (to use default values) to set
#' `hypers` argument of `calibrate()` function.
#'
#' * `dist`:
#'                       * `p1`:
#'                       * `p2`:   double (vector of doubles) representing the first parameter(s) of the
#'                                 chosen distribution(s)
#'
#' @param thetaSDist     string (vector of strings) to specify the prior distribution type(s) for
#'                       simulator correlation scale parameters
#' @param thetaSP1       double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for simulator correlation scale
#' @param thetaSP2       double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for simulator correlation scale
#' @param alphaSDist     string (vector of strings) to specify the prior distribution type(s) for
#'                       simulator correlation smoothness parameters
#' @param alphaSP1       double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for simulator correlation smoothness
#' @param alphaSP2       double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for simulator correlation smoothness
#' @param thetaBDist     string (vector of strings) to specify the prior distribution type(s) for
#'                       bias-correction correlation scale parameters
#' @param thetaBP1       double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for bias-correction correlation scale
#' @param thetaBP2       double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for bias-correction correlation scale
#' @param alphaBDist     string (vector of strings) to specify the prior distribution type(s) for
#'                       bias-correction correlation smoothness parameters
#' @param alphaBP1       double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for bias-correction correlation smoothness
#' @param alphaBP2       double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for bias-correction correlation smoothness
#' @param sigma2SDist    string (vector of strings) to specify the prior distribution type(s) for
#'                       simulator marginal variance parameter
#' @param sigma2SP1      double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for simulator marginal variance
#' @param sigma2SP2      double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for simulator marginal variance
#' @param sigma2BDist    string (vector of strings) to specify the prior distribution type(s) for
#'                       bias-correction marginal variance parameter
#' @param sigma2BP1      double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for bias-correction marginal variance
#' @param sigma2BP2      double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for bias-correction marginal variance
#' @param sigma2EDist    string (vector of strings) to specify the prior distribution type(s) for
#'                       measurement variance  parameter
#' @param sigma2EP1      double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for measurement variance
#' @param sigma2EP2      double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for measurement variance
#' @param muBDist        string (vector of strings) to specify the prior distribution type(s) for
#'                       bias-correction mean parameter
#' @param muBP1          double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for bias-correction mean
#' @param muBP2          double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for bias-correction mean
#'
#' @return nested list containing the prior specification for all calibration hyperparameters
#'
#' @example man/examples/ex_set_hyperPriors.R
#' @export
set_hyperPriors <- function(thetaSDist  = "logbeta",      thetaSP1  = 5, thetaSP2  = 5,
                            alphaSDist  = "betashift",    alphaSP1  = 5,   alphaSP2  = 2,
                            thetaBDist  = "logbeta",      thetaBP1  = 5, thetaBP2  = 5,
                            alphaBDist  = "betashift",    alphaBP1  = 5,   alphaBP2  = 2,
                            sigma2SDist = "inversegamma", sigma2SP1 = 5,   sigma2SP2 = 5,
                            sigma2BDist = "inversegamma", sigma2BP1 = 12,  sigma2BP2 = 2,
                            sigma2EDist = "inversegamma", sigma2EP1 = 1,   sigma2EP2 = 0.01,
                            muBDist     = "uniform",      muBP1     = -1,  muBP2     = 1) {
  dists <- c(thetaSDist, alphaSDist, thetaBDist, alphaBDist, sigma2SDist, sigma2BDist, sigma2EDist, muBDist)
  p1s   <- c(thetaSP1, alphaSP1, thetaBP1, alphaBP1, sigma2SP1, sigma2BP1, sigma2EP1, muBP1)
  p2s   <- c(thetaSP2, alphaSP2, thetaBP2, alphaBP2, sigma2SP2, sigma2BP2, sigma2EP2, muBP2)

  if (!is.character(dists))
    stop("'dist' argument must be character string!")
  if (sum(!(dists %in% c("uniform", "normal", "normalTr", "lognormal", "gamma", "inversegamma",
                         "beta", "betashift", "logbeta", "logistic","exponential", "fixed"))) > 0)
    stop("'dist' argument must be one of 'uniform', 'normal', 'normalTr', 'lognormal', 'gamma',
          'inversegamma', 'beta', 'betashift', 'logbeta', 'logistic','exponential', 'fixed'")
  if (sum(!is.double(p1s)) + sum(!is.double(p2s)) > 0)
      stop("Both 'p1' and `p2` arguments must be double!")
  if (length(dists) != length(p1s) || length(dists) != length(p2s))
      stop("Length of 'dist', 'p1', and 'p2' arguments must be same!")

  return(list(thetaS  = list(dist = thetaSDist,  p1 = thetaSP1,  p2 = thetaSP2),
              alphaS  = list(dist = alphaSDist,  p1 = alphaSP1,  p2 = alphaSP2),
              thetaB  = list(dist = thetaBDist,  p1 = thetaBP1,  p2 = thetaBP2),
              alphaB  = list(dist = alphaBDist,  p1 = alphaBP1,  p2 = alphaBP2),
              sigma2S = list(dist = sigma2SDist, p1 = sigma2SP1, p2 = sigma2SP2),
              sigma2B = list(dist = sigma2BDist, p1 = sigma2BP1, p2 = sigma2BP2),
              sigma2E = list(dist = sigma2EDist, p1 = sigma2EP1, p2 = sigma2EP2),
              muB     = list(dist = muBDist,     p1 = muBP1,     p2 = muBP2)))
}




