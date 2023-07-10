#' Update Augmented Covariance Matrix
#'
#' In each iteration, a new parameter value is proposed which will change some components
#' (correlation matrices and other scalers) of the augmented covariance matrix. Based on
#' the index of the changed parameter, necessary components are updated and used in
#' rebuilding of augmented covariance matrix. Since, only determinant and inverse of the
#' covariance matrix is actually used, these results are stored in cache. Note that when
#' proposed parameter is muB, the components of the augmented covariance matrix are not
#' changed.
#'
#' @param phi       vector of doubles containing the most recent update of the parameters
#' @param ichanged  index of the changed parameter in this iteration
#' @noRd
update_cov <- function(phi, ichanged) {
  if (ichanged %in% cache$ikappa) {
    cache$Xk     <- Xk <-  cbind(cache$Xf, matrix(replicate(cache$n, phi[cache$ikappa]), nrow=cache$n, byrow=T))
    cache$CorKK  <- correlation(Xk, theta = phi[cache$ithetaS], alpha = phi[cache$ialphaS])
    cache$CorKS  <- correlation(Xk, cache$Xs, phi[cache$ithetaS], phi[cache$ialphaS])
    cache$CorSK  <- t(cache$CorKS)

  } else if (ichanged %in% cache$ithetaS | ichanged %in% cache$ialphaS) {
    cache$CorKK  <- correlation(cache$Xk, theta = phi[cache$ithetaS], alpha = phi[cache$ialphaS])
    cache$CorKS  <- correlation(cache$Xk, cache$Xs, phi[cache$ithetaS], phi[cache$ialphaS])
    cache$CorSK  <- t(cache$CorKS)
    cache$CorSS  <- correlation(cache$Xs, theta = phi[cache$ithetaS], alpha = phi[cache$ialphaS])

  } else if (ichanged %in% cache$ithetaB | ichanged %in% cache$ialphaB) {
    cache$CorFF   <- correlation(cache$Xf, theta = phi[cache$ithetaB], alpha = phi[cache$ialphaB])
  }

  AugCov          <- rbind(cbind(phi[cache$isigma2S]*cache$CorKK +
                                   phi[cache$isigma2B]*cache$CorFF +
                                   phi[cache$isigma2E]*cache$Inn ,
                                 phi[cache$isigma2S]*cache$CorKS),
                           cbind(phi[cache$isigma2S]*cache$CorSK,
                                 phi[cache$isigma2S]*cache$CorSS))

  CholCov         <- chol(AugCov)
  cache$InvCov    <- chol2inv(CholCov)
  cache$logDetCov <- 2 * sum(log(diag(CholCov)))
}
