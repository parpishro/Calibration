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
#' @param thetaS  prior specification of simulator correlation scale
#' @param alphaS  prior specification of simulator correlation smoothness
#' @param thetaB  prior specification of bias-correction correlation scale
#' @param alphaB  prior specification of bias-correction correlation smoothness
#' @param sigma2S prior specification of simulator marginal variance
#' @param sigma2B prior specification of bias-correction marginal variance
#' @param sigma2E prior specification of measurement error variance
#' @param muB     prior specification of bias-correction constant mean
#'
#' @return nested list containing the prior specification for all calibration hyperparameters
#' @export
#'
#' @examples
set_hyperPriors <- function(thetaSDist  = "gamma",        thetaSInit  = 0.1, thetaSP1  = 1.5, thetaSP2  = 0.5,
                            alphaSDist  = "betashift",    alphaSInit  = 1.5, alphaSP1  = 5,   alphaSP2  = 2,
                            thetaBDist  = "gamma",        thetaBInit  = 0.1, thetaBP1  = 1.5, thetaBP2  = 0.5,
                            alphaBDist  = "betashift",    alphaBInit  = 1.5, alphaBP1  = 5,   alphaBP2  = 2,
                            sigma2SDist = "inversegamma", sigma2SInit = 1 ,  sigma2SP1 = 100, sigma2SP2 = 0.1,
                            sigma2BDist = "inversegamma", sigma2BInit = 1,   sigma2BP1 = 1,   sigma2BP2 = 1,
                            sigma2EDist = "inversegamma", sigma2EInit = 1,   sigma2EP1 = 10,  sigma2EP2 = 1,
                            muBDist     = "uniform",      muBInit     = 0,   muBP1     = -1,  muBP2     = 1) {
  stopifnot(thetaSDist %in% c("uniform", "gaussian", "gamma", "beta", "lognormal", "logistic",
                                   "betashift", "exponential", "inversegamma", "jeffreys", "fixed")  &&
              is.double(thetaSInit) && is.double(thetaSP1) && is.double(thetaSP2) &&
              length(thetaSDist) == length(thetaSP1) &&
              length(thetaSDist) == length(thetaSP2) &&
              length(thetaSDist) == length(thetaSInit))
  stopifnot(alphaSDist %in% c("uniform", "gaussian", "gamma", "beta", "lognormal", "logistic",
                                   "betashift", "exponential", "inversegamma", "jeffreys", "fixed") &&
              is.double(alphaSInit) && is.double(alphaSP1) && is.double(alphaSP2) &&
              length(alphaSDist) == length(alphaSP1) &&
              length(alphaSDist) == length(alphaSP2) &&
              length(alphaSDist) == length(alphaSInit))
  stopifnot(thetaBDist %in% c("uniform", "gaussian", "gamma", "beta", "lognormal", "logistic",
                                   "betashift", "exponential", "inversegamma", "jeffreys", "fixed") &&
              is.double(thetaBInit) && is.double(thetaBP1) && is.double(thetaBP2) &&
              length(thetaBDist) == length(thetaBP1) &&
              length(thetaBDist) == length(thetaBP2) &&
              length(thetaBDist) == length(thetaBInit))
  stopifnot(alphaBDist %in% c("uniform", "gaussian", "gamma", "beta", "lognormal", "logistic",
                               "betashift", "exponential", "inversegamma", "jeffreys", "fixed") &&
              is.double(alphaBInit) && is.double(alphaBP1) && is.double(alphaBP2) &&
              length(alphaBDist) == length(alphaBP1) &&
              length(alphaBDist) == length(alphaBP2) &&
              length(alphaBDist) == length(alphaBInit))
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




