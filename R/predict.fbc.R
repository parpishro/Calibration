#' Predict based on Calibrated Model
#'
#' It uses MCMC samples of calibration parameters and predicts the response for new input
#' configurations. For every vector of joint parameter draw, a prediction is made to get a
#' vector of prediction for each new input configuration (and a matrix of predictions for
#' all new input configurations). This vector represents the estimated distribution of
#' response conditioned on parameter samples. Full Bayesian framework of this package
#' enables also uncertainty quantification.
#'
#' @param object  class fbc (output of `calibrate()` function)
#' @param newdata matrix of new field input configuration (\eqn{\bf{X}_f^*}), where each
#'                row represents a vector of field input configuration (\eqn{\bf{x}_f^*})
#'                and columns represent experimental inputs
#'
#' @return list containing two vectors:
#'            - `pred` representing the estimate (\eqn{\bf{y}_f^*})
#'            - `se` represents the uncertainty about estimates (\eqn{\bf{\sigma^2_y}})
#'
#' @export
#'
#' @example
predict.fbc <- function(object, newdata, type="Bayesian") {
  c              <- object$cache
  Phi            <- object$Phi
  np             <- nrow(Phi)
  nx             <- nrow(newdata)
  s              <- c$scale
  preds          <- double(nx)
  vars           <- double(nx)
  predMean       <- matrix(0, nrow = np, ncol = nx)
  predVar        <- matrix(0, nrow = np, ncol = nx)
  Xstar          <- matrix(scale(newdata, center=s$expMin, scale=s$expRange), ncol=c$q)
  Phi[,c$ikappa] <- matrix(scale(Phi[,c$ikappa], center=s$calMin, scale=s$calRange), ncol=c$q)
  inds           <- c(c$indices, n=c$n, m=c$m)

  compute_prediction <- function(xstar, phi, InvCov, res) {
    xkStar <- matrix(c(xstar, phi[c$ikappa]), nrow=1)
    corSN  <- correlation(c$Xs, xkStar,                theta=phi[c$ithetaS], alpha=phi[c$ialphaS])
    corKN  <- correlation(Xk,   xkStar,                theta=phi[c$ithetaS], alpha=phi[c$ialphaS])
    corFN  <- correlation(c$Xf, matrix(xstar, ncol=c$p), theta=phi[c$ithetaB], alpha=phi[c$ialphaB])
    covND  <- matrix(c(phi[c$isigma2S]*corKN + phi[c$isigma2B]*corFN, phi[c$isigma2S]*corSN), ncol=1)

    pMean <- phi[c$imuB] + (t(covND) %*% InvCov %*% res)
    pVar  <- (phi[c$isigma2S]+phi[c$isigma2B]+phi[c$isigma2E]) - (t(covND) %*% InvCov %*% covND) + ((1-sum(InvCov%*%covND))^2)/sum(InvCov)
    return(list(pMean=pMean, pVar=pVar))
  }

  if (type == "Bayesian") {
    for (i in 1:np) {
      phi    <- as.double(Phi[i, ,drop=T])
      Xk     <- cbind(c$Xf,    matrix(as.double(replicate(c$n, phi[c$ikappa])), nrow=c$n))
      InvCov <- compute_covs(phi, c$Xf, c$Xs, Xk, inds)
      res    <- matrix(c$y - phi[c$imuB], ncol = 1)
      for (j in 1:nx) {
        prediction    <- compute_prediction(Xstar[j, ,drop=F], phi, InvCov, res)
        predMean[i,j] <- prediction$pMean
        predVar[i,j]  <- prediction$pVar
      }
    }
    predMean <- (predMean*s$sdYs)+s$meanYs
    predVar  <- predVar*(s$sdYs^2)
    varMean  <- round(apply(predVar,  2, mean), 6)
    meanVar  <- round(apply(predMean, 2, var),  6)

    preds    <- round(apply(predMean, 2, mean), 3)
    vars     <- varMean + meanVar
  } else if (type == "MAP") {
    iMAP     <- which.max(object$logPost)
    phi      <- as.double(Phi[iMAP, ,drop=T])
    Xk       <- cbind(c$Xf,    matrix(as.double(replicate(c$n, phi[c$ikappa])), nrow=c$n))
    InvCov   <- compute_covs(phi, c$Xf, c$Xs, Xk, inds)
    res      <- matrix(c$y - phi[c$imuB], ncol = 1)

    for (j in 1:nx) {
      prediction <- compute_prediction(Xstar[j, ,drop=F], phi, InvCov, res)
      preds[j]   <- prediction$pMean
      vars[j]    <- prediction$pVar
    }
    preds    <- (preds*s$sdYs)+s$meanYs
    vars     <- vars*(s$sdYs^2)
  } else {
    stop("Invalid prediction type!")
  }

  return(list(pred=preds, se=round(sqrt(vars), 4)))
}


