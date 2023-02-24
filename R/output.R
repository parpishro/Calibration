output <- function() {
  Params     <- matrix(nrow = nrow(cache$Params), ncol = cache$k)
  paramNames <- character(cache$k)
  for (i in 1:cache$k) {
    if (i %in% cache$itheta) {
      paramNames[i]        <- paste0("theta", i)
      Params[, i] <- cache$Params[i]
    } else if (i %in% cache$ilambdaS) {
      paramNames[i] <- paste0("rhoS", i - cache$ilambdaS[1] + 1)
      Params[, i] <- 1 - 1/(1+exp(cache$Params[, i]))
    } else if (i %in% cache$igammaS) {
      paramNames[i] <- paste0("alphaS", i - cache$igammaS[1] + 1)
      Params[, i] <- 1 + 1/(1+exp(-cache$Params[, i]))
    } else if (i %in% cache$ilambdaB) {
      paramNames[i] <- paste0("rhoB", i - cache$ilambdaB[1] + 1)
      Params[, i] <- 1 - 1/(1+exp(cache$Params[, i]))
    } else if (i %in% cache$igammaB) {
      paramNames[i] <- paste0("alphaB", i - cache$igammaB[1] + 1)
      Params[, i] <- 1 + 1/(1+exp(-cache$Params[, i]))
    } else if (i %in% cache$isigma2S) {
      paramNames[i] <- "sigma2S"
      Params[, i] <- abs(cache$Params[, i])
    } else if (i %in% cache$isigma2B) {
      paramNames[i] <- "sigma2B"
      Params[, i] <- abs(cache$Params[, i])
    } else if (i %in% cache$isigma2E) {
      paramNames[i] <- "sigma2E"
      Params[, i] <- abs(cache$Params[, i])
    } else if (i %in% cache$imuHat) {
      paramNames[i] <- "muHat"
      Params[, i] <- cache$Params[, i]
    }
  }
  paramMean      <- apply(Params, 2, mean)
  paramVar       <- apply(Params, 2, var)
  estimates      <- data.frame(var=paramNames, est=paramMean, se=sqrt(paramVar))
  dist           <- as.data.frame(Params)
  colnames(dist) <- paramNames
  return(list(estimates     = estimates,
              distributions = dist,
              logPost       = cache$logPost,
              vars          = paramNames,
              summary       = summary(dist)))
}
