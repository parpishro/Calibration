#' Model Output
#'
#' `output` creates a `fbc` object using all information that might be useful for user.
#' The main component of the output is matrix `Phi`, which represents draws from the joint
#' posterior distribution of all calibration model parameters.
#'
#' @return a `fbc` object containing:
#'  * Phi:            A numeric matrix in which each row represents a draw from joint
#'                    posterior distribution of parameters (after thinning to remove
#'                    autocorrelation between consecutive draws) and each column
#'                    represents a parameter of the model. In essence each column
#'                    approximates the marginal distribution of that parameter.
#'  * estimates:      a dataframe containing the summary of estimates of parameters and
#'                    their uncertainty
#'  * logPost:        a numeric vector of same length as Phi rows representing the log of
#'                    posterior likelihood
#'  * priors:         prior specifications of all parameters
#'  * acceptance:     a numeric vector representing the final acceptance rate of MH MCMC
#'                    algorithm
#'  * vars:           name of all parameters (based on below notation)
#'  * cache:          an environment containing original datasets and indexing variables
#'                    that is used in `predict()` function.
#' @import stats
#' @noRd
output <- function() {
  Params        <- cache$Params
  paramMean     <- round(apply(Params, 2, mean),   4)
  paramMedian   <- round(apply(Params, 2, median), 4)
  paramMode     <- round(apply(Params, 2, pmode),  4)
  paramSd       <- round(apply(Params, 2, sd),     4)
  param50Lwr    <- round(apply(Params, 2, quantile, 0.25), 4)
  param50Upr    <- round(apply(Params, 2, quantile, 0.75), 4)
  param80Lwr    <- round(apply(Params, 2, quantile, 0.1),  4)
  param80Upr    <- round(apply(Params, 2, quantile, 0.9),  4)
  estimates     <- data.frame(mean  = paramMean,  median = paramMedian, mode = paramMode,
                              lwr50 = param50Lwr, upr50  = param50Upr,
                              lwr80 = param80Lwr, upr80  = param80Upr,  sd   = paramSd)
  Phi           <- round(as.data.frame(Params), 4)
  colnames(Phi) <- cache$priors$param
  obj           <- list(Phi        = Phi,
                        estimates  = estimates,
                        logPost    = cache$logPost,
                        priors     = cache$priors,
                        acceptance = cache$acceptance,
                        vars       = cache$priors$param,
                        data       = list(Xf = cache$Xf, Xs = cache$Xs, y = cache$y),
                        scale      = cache$scale,
                        indices    = cache$indices)
  fbcObj        <- fbc(obj)

  return(fbcObj)
}
