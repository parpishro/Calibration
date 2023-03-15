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
#'
#' @return
update_cov <- function(phi, iChanged) {

  # indices for parameters in phi
  ikappa    <- get('ikappa',   envir = cache)
  ithetaS   <- get('ithetaS',  envir = cache)
  ialphaS   <- get('ialphaS',  envir = cache)
  ithetaB   <- get('ithetaB',  envir = cache)
  ialphaB   <- get('ialphaB',  envir = cache)
  isigma2S  <- get('isigma2S', envir = cache)
  isigma2B  <- get('isigma2B', envir = cache)
  isigma2E  <- get('isigma2E', envir = cache)

  if (iChanged %in% ikappa) {
    Xf     <- cbind(cache$Xb, replicate(cache$n, phi[ikappa], cache$n))
    CorFF  <- correlation(Xf, theta = phi[ithetaS], alpha = phi[ialphaS])
    CorFS  <- correlation(Xf, cache$Xs, phi[ithetaS], phi[ialphaS])
    CorSF  <- t(CorFS)
    assign('Xf',    Xf,    envir = cache)
    assign('CorFF', CorFF, envir = cache)
    assign('CorFS', CorFS, envir = cache)
    assign('CorSF', CorSF, envir = cache)

  } else if (iChanged %in% ithetaS | iChanged %in% ialphaS) {
    CorFF  <- correlation(cache$Xf, theta = phi[ithetaS], alpha = phi[ialphaS])
    CorFS  <- correlation(cache$Xf, cache$Xs, phi[ithetaS], phi[ialphaS])
    CorSF  <- t(CorFS)
    CorSS  <- correlation(cache$Xs, theta = phi[ithetaS], alpha = phi[ialphaS])
    assign('CorSS', CorSS, envir = cache)
    assign('CorFF', CorFF, envir = cache)
    assign('CorFS', CorFS, envir = cache)
    assign('CorSF', CorSF, envir = cache)

  } else if (iChanged %in% ithetaB | iChanged %in% ialphaB) {
    CorB   <- correlation(cache$Xb, theta = phi[ithetaB], alpha = phi[ialphaB])
    assign('CorB', CorB, envir = cache)

  } else if (iChanged %in% isigma2S)
    assign('sigma2S', phi[isigma2S], envir = cache)

  else if (iChanged %in% isigma2B)
    assign('sigma2B', phi[isigma2B], envir = cache)

  else if (iChanged %in% isigma2E)
    assign('sigma2E', phi[isigma2E], envir = cache)

  else stop("invalid iChanged argument!")

  Inn     <- diag(cache$n)
  AugCov  <- rbind(cbind(cache$sigma2S*cache$CorFF + cache$sigma2B*cache$CorB + cache$sigma2E*Inn, cache$sigma2S*cache$CorFS),
                   cbind(cache$sigma2S*cache$CorSF, cache$sigma2S*cache$CorSS))
  CholCov <- try(chol(AugCov), silent = T)

  if (length(class(CholCov)) == 1 && class(CholCov) == "try-error")
    return(NULL)

  InvCov    <- chol2inv(CholCov)
  logDetCov <- sum(2*log(diag(CholCov)))
  muHat     <- update_mu()
  res       <- cache$y - muHat
  assign('InvCov',    InvCov,    envir = cache)
  assign('logDetCov', logDetCov, envir = cache)
  assign('res',       res,       envir = cache)
  assign('muHat',     muHat,     envir = cache)


  return(list(InvCov = InvCov, logDetCov = logDetCov))

}
