output <- function() {
  Params     <- matrix(nrow = nrow(cache$Params), ncol = cache$k)
  paramNames <- character(cache$k)
  dist       <- as.data.frame(cache$Params)
  for (i in 1:cache$k) {
    if (i %in% cache$ikappa) {
      paramNames[i] <- paste0("kappa", i)
      dist[, i]     <- (dist[, i] * (cache$calibMax[i] - cache$calibMin[i])) + cache$calibMin[i]
    } else if (i %in% cache$ithetaS) {
      paramNames[i] <- paste0("thetaS", i - cache$ithetaS[1] + 1)
    } else if (i %in% cache$ialphaS) {
      paramNames[i] <- paste0("alphaS", i - cache$ialphaS[1] + 1)
    } else if (i %in% cache$ithetaB) {
      paramNames[i] <- paste0("thetaB", i - cache$ithetaB[1] + 1)
    } else if (i %in% cache$ialphaB) {
      paramNames[i] <- paste0("alphaB", i - cache$ialphaB[1] + 1)
    } else if (i %in% cache$isigma2S) {
      paramNames[i] <- "sigma2S"
    } else if (i %in% cache$isigma2B) {
      paramNames[i] <- "sigma2B"
    } else if (i %in% cache$isigma2E) {
      paramNames[i] <- "sigma2E"
    } else if (i %in% cache$imuHat) {
      paramNames[i] <- "muHat"
    }
  }

  paramMean        <- round(apply(dist, 2, mean), 2)
  paramVar         <- round(apply(dist, 2, var), 2)
  estimates        <- data.frame(var=paramNames, est=paramMean, se=paramVar)
  colnames(dist)   <- paramNames
  return(list(estimates     = estimates,
              distributions = round(dist, 2),
              logPost       = cache$logPost,
              vars          = paramNames,
              summary       = summary(dist)))
}
