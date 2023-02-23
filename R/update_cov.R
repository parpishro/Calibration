update_cov <- function(phi, changed) {
  # parameter indices in phi vector

  # indices for parameters in phi
  # indices for parameters in phi
  itheta    <- get('itheta',   envir = cache)
  ilambdaS   <- get('ilambdaS',  envir = cache)
  igammaS   <- get('igammaS',  envir = cache)
  ilambdaB   <- get('ilambdaB',  envir = cache)
  igammaB   <- get('igammaB',  envir = cache)
  isigma2S  <- get('isigma2S', envir = cache)
  isigma2B  <- get('isigma2B', envir = cache)
  isigma2E  <- get('isigma2E', envir = cache)

  # corFF : (n * n) correlation matrix of augmented Xf's
  # corFS : (n * m) correlation matrix between Xf's, Xs's
  # corSF : (m * n) correlation matrix between Xs's, Xf's
  # corSS : (m * m) correlation matrix between Xs's
  # corB  : (n * n) correlation matrix between Xb's

  if (changed %in% itheta) {
    Xf     <- cbind(cache$Xb, replicate(n, phi[itheta], n))
    CorFF  <- correlation(cache$Xf, lambda = phi[ilambdaS], gamma = phi[igammaS])
    CorFS  <- correlation(cache$Xf, cache$Xs, phi[ilambdaS], phi[igammaS])
    CorSF  <- t(CorFS)

    assign('Xf',    Xf,    envir = cache)
    assign('CorFF', CorFF, envir = cache)
    assign('CorFS', CorFS, envir = cache)
    assign('CorSF', CorSF, envir = cache)

  } else if ((changed %in% ilambdaS) || (changed %in% igammaS)) {
    CorFF  <- correlation(cache$Xf, lambda = phi[ilambdaS], gamma = phi[igammaS])
    CorFS  <- correlation(cache$Xf, cache$Xs, phi[ilambdaS], phi[igammaS])
    CorSF  <- t(CorFS)
    CorSS  <- correlation(cache$Xs, lambda = phi[ilambdaS], gamma = phi[igammaS])
    muHat  <- mu_hat()

    assign('CorSS', CorSS, envir = cache)
    assign('CorFF', CorFF, envir = cache)
    assign('CorFS', CorFS, envir = cache)
    assign('CorSF', CorSF, envir = cache)
    assign('muHat', muHat, envir = cache)

  } else if ((changed %in% ilambdaB) || (changed %in% igammaB)) {
    CorB   <- correlation(cache$Xb, lambda = phi[ilambdaB], gamma = phi[igammaB])
    assign('CorB', CorB, envir = cache)

  } else if (changed %in% isigma2S) {
    assign('sigma2S', phi[isigma2S], envir = cache)

  } else if (changed %in% isigma2B) {
    assign('sigma2B', phi[isigma2B], envir = cache)

  } else if (changed %in% isigma2E) {
    assign('sigma2E', phi[isigma2E], envir = cache)


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

  assign('Chol', CholCov, envir = cache)

  return(list(InvCov = InvCov, logDetCov = logDetCov))

}
