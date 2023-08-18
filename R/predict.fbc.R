#' Predict based on Calibrated Model
#'
#' It uses MCMC samples of calibration parameters and predicts the response for new input
#' configurations. For every vector of joint parameter draw, a prediction is made to get a
#' vector of prediction for each new input configuration (and a matrix of predictions for
#' all new input configurations). This vector represents the estimated distribution of
#' response conditioned on parameter samples. Full Bayesian framework of this package
#' enables also uncertainty quantification.
#'
#' @param object     class fbc (output of `calibrate()` function)
#' @param newdata matrix of new field input configuration (\eqn{\bf{X}_f^*}), where each
#'                row represents a vector of field input configuration (\eqn{\bf{x}_f^*})
#'                and columns represent experimental inputs
#' @param method character string from "MAP" or "Bayesian" (default) as method of prediction
#' @param ...     other arguments
#'
#' @return list containing two vectors:
#'            - `pred` representing the estimate (\eqn{\bf{y}_f^*})
#'            - `se` represents the uncertainty about estimates (\eqn{\bf{\sigma^2_y}})
#' @importFrom stats var
#' @example man/examples/ex_predict.R
#' @export
predict.fbc <- function(object, newdata, method = "Bayesian", ...) {
  Phi                 <- object$Phi
  indices             <- object$indices
  scale               <- object$scale

  np                  <- nrow(Phi)
  nx                  <- nrow(newdata)
  preds               <- double(nx)
  vars                <- double(nx)
  predMean            <- matrix(0, nrow = np, ncol = nx)
  predVar             <- matrix(0, nrow = np, ncol = nx)
  Xstar               <- matrix(scale(newdata,
                                      center = scale$expMin,
                                      scale  = scale$expRange),
                                ncol = indices$q)
  Phi[,indices$ikappa] <- matrix(scale(Phi[,indices$ikappa],
                                      center = scale$calMin,
                                      scale  = scale$calRange),
                                ncol = indices$q)

  compute_prediction <- function(xstar, phi, InvCov, res) {
    xkStar <- matrix(c(xstar, phi[indices$ikappa]), nrow = 1)
    corSN  <- correlation(object$data$Xs, xkStar,
                          theta = phi[indices$ithetaS], alpha = phi[indices$ialphaS])
    corKN  <- correlation(Xk, xkStar,
                          theta = phi[indices$ithetaS], alpha = phi[indices$ialphaS])
    corFN  <- correlation(object$data$Xf, matrix(xstar, ncol = indices$p),
                          theta = phi[indices$ithetaB], alpha = phi[indices$ialphaB])
    covND  <- matrix(c(phi[indices$isigma2S]*corKN + phi[indices$isigma2B]*corFN,
                       phi[indices$isigma2S]*corSN), ncol = 1)

    pMean <- phi[indices$imuB] + (t(covND) %*% InvCov %*% res)
    pVar  <- (phi[indices$isigma2S] + phi[indices$isigma2B] + phi[indices$isigma2E]) -
               (t(covND) %*% InvCov %*% covND) +
               ((1 - sum(InvCov %*% covND))^2)/sum(InvCov)
    return(list(pMean = pMean, pVar = pVar))
  }

  if (method == "Bayesian") {
    for (i in 1:np) {
      phi    <- as.double(Phi[i, ,drop = T])
      Xk     <- cbind(object$data$Xf,
                      matrix(as.double(replicate(indices$n, phi[indices$ikappa])), nrow = indices$n))
      InvCov <- compute_covs(phi, object$data$Xf, object$data$Xs, Xk, indices)
      res    <- matrix(object$data$y - phi[indices$imuB], ncol = 1)
      for (j in 1:nx) {
        prediction    <- compute_prediction(Xstar[j, ,drop = F], phi, InvCov, res)
        predMean[i,j] <- prediction$pMean
        predVar[i,j]  <- prediction$pVar
      }
    }
    predMean <- (predMean*scale$sdYs) + scale$meanYs
    predVar  <- predVar*(scale$sdYs^2)
    varMean  <- round(apply(predVar,  2, mean), 6)
    meanVar  <- round(apply(predMean, 2, var),  6)

    preds    <- round(apply(predMean, 2, mean), 3)
    vars     <- varMean + meanVar
  } else if (method == "MAP") {
    iMAP     <- which.max(object$logPost)
    phi      <- as.double(Phi[iMAP, ,drop = T])
    Xk       <- cbind(object$data$Xf,
                      matrix(as.double(replicate(indices$n, phi[indices$ikappa])), nrow = indices$n))
    InvCov   <- compute_covs(phi, object$data$Xf, object$data$Xs, Xk, indices)
    res      <- matrix(object$data$y - phi[indices$imuB], ncol = 1)

    for (j in 1:nx) {
      prediction <- compute_prediction(Xstar[j, ,drop = F], phi, InvCov, res)
      preds[j]   <- prediction$pMean
      vars[j]    <- prediction$pVar
    }
    preds    <- (preds*scale$sdYs) + scale$meanYs
    vars     <- vars*(scale$sdYs^2)
  } else {
    stop("Invalid prediction method!")
  }

  return(list(pred = preds, se = round(sqrt(vars), 4)))
}


