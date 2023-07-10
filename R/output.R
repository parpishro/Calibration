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
  paramNames   <- character(cache$l)
  Params       <- as.data.frame(cache$Params)
  scales       <- cache$scale
  for (i in 1:(cache$l-4)) {
    if (i %in% cache$ikappa) {
      paramNames[i] <- paste0("kappa", i)
      Params[, i]   <- Params[, i]*scales$calRange[i] + scales$calMin[i]
    } else if (i %in% cache$ithetaS) {
      paramNames[i] <- paste0("thetaS", i - cache$ithetaS[1] + 1)
    } else if (i %in% cache$ialphaS) {
      paramNames[i] <- paste0("alphaS", i - cache$ialphaS[1] + 1)
    } else if (i %in% cache$ithetaB) {
      paramNames[i] <- paste0("thetaB", i - cache$ithetaB[1] + 1)
    } else if (i %in% cache$ialphaB) {
      paramNames[i] <- paste0("alphaB", i - cache$ialphaB[1] + 1)
    }
  }

  paramNames[(cache$l-3):cache$l] <- c("sigma2S", "sigma2B", "sigma2E", "muB")
  colnames(Params) <- paramNames

  paramMean        <- round(apply(Params, 2, mean),   3)
  paramMedian      <- round(apply(Params, 2, median), 3)
  paramMode        <- round(apply(Params, 2, pmode),  3)
  paramSd          <- round(apply(Params, 2, sd),     3)
  param50Lwr       <- round(apply(Params, 2, quantile, 0.25), 3)
  param50Upr       <- round(apply(Params, 2, quantile, 0.75), 3)
  param80Lwr       <- round(apply(Params, 2, quantile, 0.1),  3)
  param80Upr       <- round(apply(Params, 2, quantile, 0.9),  3)
  estimates        <- data.frame(mean=paramMean, median = paramMedian, mode = paramMode,
                                 lwr50 = param50Lwr, upr50 = param50Upr,
                                 lwr80 = param80Lwr, upr80 = param80Upr, sd = paramSd)
  obj  <- list(Phi        = round(Params, 2),
               estimates  = estimates,
               logPost    = cache$logPost,
               priors     = cache$priors,
               acceptance = cache$acceptance,
               vars       = paramNames,
               data       = list(Xf = cache$Xf, Xs = cache$Xs, y = cache$y),
               scale      = cache$scale,
               indices    = cache$indices,
               priorFns   = cache$priorFns,
               proposalSD = cache$sdRates)
  fbcObj <- fbc(obj)

  return(fbcObj)
}
