update_cov <- function(phi, changed) {
  # parameter indices in phi vector

  # indices for parameters in phi
  # indices for parameters in phi
  iTheta    <- get('iTheta',   envir = cache)
  iOmegaS   <- get('iOmegaS',  envir = cache)
  iAlphaS   <- get('iAlphaS',  envir = cache)
  iOmegaB   <- get('iOmegaB',  envir = cache)
  iAlphaB   <- get('iAlphaB',  envir = cache)
  iSigma2S  <- get('iSigma2S', envir = cache)
  iSigma2B  <- get('iSigma2B', envir = cache)
  iSigma2E  <- get('iSigma2E', envir = cache)

  # corFF : (n * n) correlation matrix of augmented Xf's
  # corFS : (n * m) correlation matrix between Xf's, Xs's
  # corSF : (m * n) correlation matrix between Xs's, Xf's
  # corSS : (m * m) correlation matrix between Xs's
  # corB  : (n * n) correlation matrix between Xb's

  if (changed %in% iTheta) {
    Xf     <- cbind(cache$Xb, replicate(n, phi[iTheta], n))
    CorFF  <- corelation(cache$Xf, scale = phi[iOmegaS], smooth = phi[iAlphaS])
    CorFS  <- corelation(cache$Xf, cache$Xs, phi[iOmegaS], phi[iAlphaS])
    CorSF  <- t(CorFS)

    assign('Xf',    Xf,    envir = cache)
    assign('CorFF', CorFF, envir = cache)
    assign('CorFS', CorFS, envir = cache)
    assign('CorSF', CorSF, envir = cache)

  } else if ((changed %in% iOmegaS) || (changed %in% iAlphaS)) {
    CorFF  <- corelation(cache$Xf, scale = phi[iOmegaS], smooth = phi[iAlphaS])
    CorFS  <- corelation(cache$Xf, cache$Xs, phi[iOmegaS], phi[iAlphaS])
    CorSF  <- t(CorFS)
    CorSS  <- corelation(cache$Xs, scale = phi[iOmegaS], smooth = phi[iAlphaS])
    muHat  <- mu_hat(CorSS, cache$ys)
    res    <- cache$y - muHat

    assign('CorSS', CorSS, envir = cache)
    assign('CorFF', CorFF, envir = cache)
    assign('CorFS', CorFS, envir = cache)
    assign('CorSF', CorSF, envir = cache)
    assign('muHat', muHat, envir = cache)
    assign('res',   res,   envir = cache)

  } else if ((changed %in% iOmegaB) || (changed %in% iAlphaB)) {
    CorB   <- corelation(cache$Xb, scale = phi[iOmegaB], smooth = phi[iAlphaB])
    assign('CorB', CorB, envir = cache)

  } else if (changed %in% iSigma2S) {
    assign('sigma2S', phi[iSigma2S], envir = cache)

  } else if (changed %in% iSigma2B) {
    assign('sigma2B', phi[iSigma2B], envir = cache)

  } else if (changed %in% iSigma2E) {
    assign('sigma2E', phi[iSigma2E], envir = cache)


  } else stop("invalid changed argument!")

  Inn       <- diag(n)
  AugCov    <- rbind(cbind(((cache$sigma2S * cache$CorFF) +
                            (cache$sigma2B * cache$CorB) +
                            (cache$sigma2E * Inn)),
                           (cache$sigma2S * cache$CorFS)),
                     cbind((cache$sigma2S * cache$CorSF),
                           (cache$sigma2S * cache$CorSS)))
  CholCov   <- chol(AugCov)
  InvCov    <- chol2inv(CholCov)
  logDetCov <- sum(2*log(diag(CholCov)))

  return(list(InvCov = InvCov, logDetCov = logDetCov))

}
