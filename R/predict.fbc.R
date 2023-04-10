#' predict
#'
#' @param object     class fbc
#' @param newdata a vector of new data
#'
#' @return predictions
#'
#' @export
predict.fbc <- function(object, newdata) {
  cache <- object$cache
  phi   <- object$estimates$mean
  cache$Xf     <- Xf <-  cbind(cache$Xb, replicate(cache$n, phi[cache$ikappa], cache$n))
  cache$CorFF  <- correlation(Xf, theta = phi[cache$ithetaS], alpha = phi[cache$ialphaS])
  cache$CorFS  <- correlation(Xf, cache$Xs, phi[cache$ithetaS], phi[cache$ialphaS])
  cache$CorSF  <- t(cache$CorFS)
  cache$CorSS  <- correlation(cache$Xs, theta = phi[cache$ithetaS], alpha = phi[cache$ialphaS])
  cache$CorB   <- correlation(cache$Xb, theta = phi[cache$ithetaB], alpha = phi[cache$ialphaB])

  Inn     <- diag(cache$n)
  AugCov  <- rbind(cbind(phi[cache$isigma2S]*cache$CorFF + phi[cache$isigma2B]*cache$CorB + phi[cache$isigma2E]*Inn, phi[cache$isigma2S]*cache$CorFS),
                   cbind(phi[cache$isigma2S]*cache$CorSF, phi[cache$isigma2S]*cache$CorSS))
  CholCov <- try(chol(AugCov), silent = T)

  if (length(class(CholCov)) == 1 && class(CholCov) == "try-error")
    return(NULL)

  InvCov    <- chol2inv(CholCov)
  logDetCov <- 2 * sum(log(diag(CholCov)))


  x                                    <- scale(newdata, center=cache$experMin, scale=cache$experMax-cache$experMin)
  xstar                                <- double(cache$p+cache$q)
  xstar[1:cache$p]                     <- scale(newdata, center=cache$experMin, scale=cache$experMax-cache$experMin)
  xstar[(cache$p+1):(cache$p+cache$q)] <- scale(phi[cache$ikappa], center=cache$calibMin, scale=cache$calibMax-cache$calibMin)



  covVec   <- phi[cache$isigma2S]*rbind(correlation(cache$Xf, matrix(xstar, nrow = 1), theta = phi[cache$ithetaS], alpha = phi[cache$ialphaS]),
                                        correlation(cache$Xs, matrix(xstar, nrow = 1), theta = phi[cache$ithetaS], alpha = phi[cache$ialphaS])) +
              phi[cache$isigma2B]*rbind(correlation(cache$Xb, matrix(x, nrow = 1), theta = phi[cache$ithetaB], alpha = phi[cache$ialphaB]),
                                        matrix(0, nrow = cache$m, ncol = 1))

  preds    <- t(covVec) %*% InvCov %*% cache$y
  return((preds* sd(cache$ys)) + mean(cache$ys))
}
fbc <- setClass("fbc")
setMethod("predict", signature = "fbc", definition = predict.fbc)
