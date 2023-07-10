#' Set Calibration Hyperparameter Priors
#'
#' `set_hyperPriors()` can be used to specify one or more hyperparameter priors. There are
#' 8 classes of hyperparameters. Each class can be specified using a named list with four
#' fields:
#'  - `dist`: character string (vector of character strings) to specify the family of the
#'            distributions. All common probability distribution are implemented.
#'  - `init`: initial value to start the MCMC algorithm
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
#'                       * `init`:
#'
#' @param thetaSDist     string (vector of strings) to specify the prior distribution type(s) for
#'                       simulator correlation scale parameters
#' @param thetaSInit     double (vector of doubles) that represent initial value(s) of simulator
#'                       correlation scale parameter parameters to start the MCMC
#' @param thetaSP1       double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for simulator correlation scale
#' @param thetaSP2       double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for simulator correlation scale
#' @param alphaSDist     string (vector of strings) to specify the prior distribution type(s) for
#'                       simulator correlation smoothness parameters
#' @param alphaSInit     double (vector of doubles) that represent initial value(s) of simulator
#'                       correlation smoothness parameter parameters to start the MCMC
#' @param alphaSP1       double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for simulator correlation smoothness
#' @param alphaSP2       double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for simulator correlation smoothness
#' @param thetaBDist     string (vector of strings) to specify the prior distribution type(s) for
#'                       bias-correction correlation scale parameters
#' @param thetaBInit     double (vector of doubles) that represent initial value(s) of
#'                       bias-correction correlation scale parameter parameters to start the MCMC
#' @param thetaBP1       double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for bias-correction correlation scale
#' @param thetaBP2       double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for bias-correction correlation scale
#' @param alphaBDist     string (vector of strings) to specify the prior distribution type(s) for
#'                       bias-correction correlation smoothness parameters
#' @param alphaBInit     double (vector of doubles) that represent initial value(s) of
#'                       bias-correction  correlation smoothness parameter parameters to start the MCMC
#' @param alphaBP1       double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for bias-correction correlation smoothness
#' @param alphaBP2       double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for bias-correction correlation smoothness
#' @param sigma2SDist    string (vector of strings) to specify the prior distribution type(s) for
#'                       simulator marginal variance parameter
#' @param sigma2SInit    double (vector of doubles) that represent initial value(s) of simulator
#'                       marginal variance parameter parameters to start the MCMC
#' @param sigma2SP1      double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for simulator marginal variance
#' @param sigma2SP2      double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for simulator marginal variance
#' @param sigma2BDist    string (vector of strings) to specify the prior distribution type(s) for
#'                       bias-correction marginal variance parameter
#' @param sigma2BInit    double (vector of doubles) that represent initial value(s) of
#'                       bias-correction marginal variance parameter parameters to start the MCMC
#' @param sigma2BP1      double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for bias-correction marginal variance
#' @param sigma2BP2      double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for bias-correction marginal variance
#' @param sigma2EDist    string (vector of strings) to specify the prior distribution type(s) for
#'                       measurement variance  parameter
#' @param sigma2EInit    double (vector of doubles) that represent initial value(s) of measurement
#'                       variance parameter parameters to start the MCMC
#' @param sigma2EP1      double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for measurement variance
#' @param sigma2EP2      double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for measurement variance
#' @param muBDist        string (vector of strings) to specify the prior distribution type(s) for
#'                       bias-correction mean parameter
#' @param muBInit        double (vector of doubles) that represent initial value(s) of
#'                       bias-correction mean parameter parameters to start the MCMC
#' @param muBP1          double (vector of doubles) representing the first parameter(s) of the
#'                       chosen distribution(s) for bias-correction mean
#' @param muBP2          double (vector of doubles) representing the second parameter(s) of the
#'                       chosen distribution(s) for bias-correction mean
#'
#' @return nested list containing the prior specification for all calibration hyperparameters
#'
#' @example man/examples/ex_set_hyperPriors.R
#' @export
set_hyperPriors <- function(thetaSDist  = "gamma",        thetaSInit  = 0.1, thetaSP1  = 1.5, thetaSP2  = 0.5,
                            alphaSDist  = "betashift",    alphaSInit  = 1.5, alphaSP1  = 5,   alphaSP2  = 2,
                            thetaBDist  = "gamma",        thetaBInit  = 0.1, thetaBP1  = 1.5, thetaBP2  = 0.5,
                            alphaBDist  = "betashift",    alphaBInit  = 1.5, alphaBP1  = 5,   alphaBP2  = 2,
                            sigma2SDist = "inversegamma", sigma2SInit = 1 ,  sigma2SP1 = 100, sigma2SP2 = 0.1,
                            sigma2BDist = "inversegamma", sigma2BInit = 1,   sigma2BP1 = 1,   sigma2BP2 = 1,
                            sigma2EDist = "inversegamma", sigma2EInit = 1,   sigma2EP1 = 10,  sigma2EP2 = 1,
                            muBDist     = "uniform",      muBInit     = 0,   muBP1     = -1,  muBP2     = 1) {
  stopifnot(sum(!(thetaSDist %in% c("uniform", "gaussian", "gamma", "beta", "lognormal", "logistic",
                                    "betashift", "exponential", "inversegamma", "jeffreys", "fixed"))) == 0  &&
              is.double(thetaSInit) && is.double(thetaSP1) && is.double(thetaSP2) &&
              length(thetaSDist) == length(thetaSP1) &&
              length(thetaSDist) == length(thetaSP2))
  stopifnot(sum(!(alphaSDist %in% c("uniform", "gaussian", "gamma", "beta", "lognormal", "logistic",
                                    "betashift", "exponential", "inversegamma", "jeffreys", "fixed"))) == 0 &&
              is.double(alphaSInit) && is.double(alphaSP1) && is.double(alphaSP2) &&
              length(alphaSDist) == length(alphaSP1) &&
              length(alphaSDist) == length(alphaSP2))
  stopifnot(sum(!(thetaBDist %in% c("uniform", "gaussian", "gamma", "beta", "lognormal", "logistic",
                                    "betashift", "exponential", "inversegamma", "jeffreys", "fixed"))) == 0 &&
              is.double(thetaBInit) && is.double(thetaBP1) && is.double(thetaBP2) &&
              length(thetaBDist) == length(thetaBP1) &&
              length(thetaBDist) == length(thetaBP2))
  stopifnot(sum(!(alphaBDist %in% c("uniform", "gaussian", "gamma", "beta", "lognormal", "logistic",
                                    "betashift", "exponential", "inversegamma", "jeffreys", "fixed"))) == 0 &&
              is.double(alphaBInit) && is.double(alphaBP1) && is.double(alphaBP2) &&
              length(alphaBDist) == length(alphaBP1) &&
              length(alphaBDist) == length(alphaBP2))
  stopifnot(sigma2SDist %in% c("uniform", "gaussian", "gamma", "beta", "lognormal", "logistic",
                                   "betashift", "exponential", "inversegamma", "jeffreys", "fixed")  &&
              is.double(sigma2SInit) && is.double(sigma2SP1) && is.double(sigma2SP2))
  stopifnot(sigma2BDist %in% c("uniform", "gaussian", "gamma", "beta", "lognormal", "logistic",
                                   "betashift", "exponential", "inversegamma", "jeffreys")  &&
              is.double(sigma2BInit) && is.double(sigma2BP1) && is.double(sigma2BP2))
  stopifnot(sigma2EDist %in% c("uniform", "gaussian", "gamma", "beta", "lognormal", "logistic",
                                   "betashift", "exponential", "inversegamma", "jeffreys", "fixed")  &&
              is.double(sigma2EInit) && is.double(sigma2EP1) && is.double(sigma2EP2))
  stopifnot(muBDist %in% c("uniform", "gaussian", "gamma", "beta", "lognormal", "logistic",
                                   "betashift", "exponential", "inversegamma", "jeffreys", "fixed")  &&
              is.double(muBInit) && is.double(muBP1) && is.double(muBP2))
  return(list(thetaS  = list(dist = thetaSDist,  init = thetaSInit,  p1 = thetaSP1,  p2 = thetaSP2),
              alphaS  = list(dist = alphaSDist,  init = alphaSInit,  p1 = alphaSP1,  p2 = alphaSP2),
              thetaB  = list(dist = thetaBDist,  init = thetaBInit,  p1 = thetaBP1,  p2 = thetaBP2),
              alphaB  = list(dist = alphaBDist,  init = alphaBInit,  p1 = alphaBP1,  p2 = alphaBP2),
              sigma2S = list(dist = sigma2SDist, init = sigma2SInit, p1 = sigma2SP1, p2 = sigma2SP2),
              sigma2B = list(dist = sigma2BDist, init = sigma2BInit, p1 = sigma2BP1, p2 = sigma2BP2),
              sigma2E = list(dist = sigma2EDist, init = sigma2EInit, p1 = sigma2EP1, p2 = sigma2EP2),
              muB     = list(dist = muBDist,     init = muBInit,     p1 = muBP1,     p2 = muBP2)))
}




