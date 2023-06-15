#' Model Output
#'
#' `output` loads all of the information that might be useful for user into a
#' single list, which will be returned to user by `calibrate`
#'
#' @return a list containing the parameters estimates based on MCMCM sampling results,
#' their distribution, acceptence rate of the algorithm, log posterior values of the
#' MCMC runs, raw unfiltered result of sampling (without burn-in or thining), parameter
#' names and a quantile summary of each parameter
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

  paramNames[(cache$l-3):cache$l] <- c("sigma2S", "sigma2B", "sigma2E", "mu")
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
  obj  <- new("fbc",
              estimates  = estimates,
              Phi        = round(Params, 2),
              logPost    = cache$logPost,
              acceptance = cache$acceptance,
              vars       = paramNames,
              cache      = cache)
  return(obj)
}

pmode <- function(vec) {
  bins   <- seq(min(vec), max(vec), length.out = 101)
  counts <- double(100)
  for (i in 2:101)
    counts[i-1] <- sum(vec >= bins[i-1] & vec < bins[i])
  return(mean(bins[which.max(counts):(which.max(counts)+1)]))
}
