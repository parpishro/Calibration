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
set_hyperPriors <- function(thetaS  = list(dist="gamma",        init=0.5, p1=1.1, p2=0.1),
                            alphaS  = list(dist="betashift",    init=1.7, p1=5,   p2=2),
                            thetaB  = list(dist="gamma",        init=0.5, p1=1.1, p2=0.1),
                            alphaB  = list(dist="betashift",    init=1.7, p1=5,   p2=2),
                            sigma2S = list(dist="inversegamma", init=1 ,  p1=0.01, p2=0.01),
                            sigma2B = list(dist="inversegamma", init=1,   p1=0.01, p2=0.01),
                            sigma2E = list(dist="inversegamma", init=1,   p1=0.01, p2=0.01),
                            muB     = list(dist="uniform",      init=0,   p1=-5,  p2=5)) {
  return(list(thetaS  = thetaS,  alphaS  = alphaS,
              thetaB  = thetaB,  alphaB  = alphaB,
              sigma2S = sigma2S, sigma2B = sigma2B,
              sigma2E = sigma2E, muB     = muB))
}


