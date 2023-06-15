#' Update the augmented covariance matrix in each iteration of MCMC
#'
#' In each iteration, a new parameter value is proposed which will change some
#' components (correlation matrices and other scalers) of the augmented covariance
#' matrix. Based on the category of the changed parameter, the corresponding
#' components are updated which requires rebuild of augmented covariance matrix.
#' In the algorithm, determinant and inverse of the covariance matrix is actually used.
#' Therefore, only determinant and inverse are returned and updated components are
#' overwrite their previous value in the cache.
#'
#' @param phi      a vector of doubles containing the most recent update of the parameters
#' @param iChanged      index of the changed parameter in this iteration
#' @noRd
#' @return A list consisting of the inverse covariance matrix and the determinant of covariance matrix
update_cov <- function(phi, iChanged) {
  if (iChanged %in% cache$ikappa) {
    cache$Xk     <- Xk <-  cbind(cache$Xf, matrix(replicate(cache$n, phi[cache$ikappa]), nrow=cache$n, byrow=T))
    cache$CorKK  <- correlation(Xk, theta = phi[cache$ithetaS], alpha = phi[cache$ialphaS])
    cache$CorKS  <- correlation(Xk, cache$Xs, phi[cache$ithetaS], phi[cache$ialphaS])
    cache$CorSK  <- t(cache$CorKS)

  } else if (iChanged %in% cache$ithetaS | iChanged %in% cache$ialphaS) {
    cache$CorKK  <- correlation(cache$Xk, theta = phi[cache$ithetaS], alpha = phi[cache$ialphaS])
    cache$CorKS  <- correlation(cache$Xk, cache$Xs, phi[cache$ithetaS], phi[cache$ialphaS])
    cache$CorSK  <- t(cache$CorKS)
    cache$CorSS  <- correlation(cache$Xs, theta = phi[cache$ithetaS], alpha = phi[cache$ialphaS])

  } else if (iChanged %in% cache$ithetaB | iChanged %in% cache$ialphaB) {
    cache$CorFF   <- correlation(cache$Xf, theta = phi[cache$ithetaB], alpha = phi[cache$ialphaB])

  }



  AugCov  <- rbind(cbind(phi[cache$isigma2S]*cache$CorKK +
                         phi[cache$isigma2B]*cache$CorFF +
                         phi[cache$isigma2E]*cache$Inn ,
                         phi[cache$isigma2S]*cache$CorKS),
                   cbind(phi[cache$isigma2S]*cache$CorSK,
                         phi[cache$isigma2S]*cache$CorSS))

  CholCov <- try(chol(AugCov), silent = T)

  if (length(class(CholCov)) == 1 && class(CholCov) == "try-error")
    return(NULL)

  InvCov    <- chol2inv(CholCov)
  logDetCov <- 2 * sum(log(diag(CholCov)))

  return(list(InvCov = InvCov, logDetCov = logDetCov))

}
