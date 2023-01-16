update_cov <- function(phi, changed, env) {

  env0             <- environment()
  parent.env(env0) <- env

  # parameter indices in phi vector



  # corFF : (n * n) correlation matrix of augmented Xf's
  # corFS : (n * m) correlation matrix between Xf's, Xs's
  # corSF : (m * n) correlation matrix between Xs's, Xf's
  # corSS : (m * m) correlation matrix between Xs's
  # corB  : (n * n) correlation matrix between Xb's

  if (changed %in% calib) {
    Xf     <- cbind(Xb, replicate(n, phi[calib], n))
    CorFF  <- corelation(Xf, Xf, scale = phi[scaleS], smooth = phi[smoothS])
    CorFS  <- corelation(Xf, Xs, scale = phi[scaleS], smooth = phi[smoothS])
    CorSF  <- t(CorFS)

  } else if ((changed %in% scaleS) || (changed %in% smoothS)) {
    CorFF  <- corelation(Xf, Xf, scale = phi[scaleS], smooth = phi[smoothS])
    CorFS  <- corelation(Xf, Xs, scale = phi[scaleS], smooth = phi[smoothS])
    CorSF  <- t(CorFS)
    CorSS  <- corelation(Xs, Xs, scale = phi[scaleS], smooth = phi[smoothS])
    muHat  <- mu_hat(CorSS, ys)
    res    <- y - muHat

  } else if ((changed %in% scaleB) || (changed %in% smoothB)) {
    CorB   <- corelation(Xb, Xb, scale = phi[scaleB], smooth = phi[smoothB])

  } else if (changed %in% sig2S) {
    sig2S  <- phi[sig2S]

  } else if (changed %in% sig2B) {
    sig2B  <- phi[sig2B]

  } else if (changed %in% sig2E) {
    sig2E  <- phi[sig2E]

  } else stop("invalid changed argument!")

  I        <- diag(n)
  AugCov   <- rbind(cbind(((sig2S * CorFF) + (sig2B * CorB) + (sig2E * I)), (sig2S * CorFS)),
                    cbind((sig2S * CorSF), (sig2S * CorSS)))

  return(chol_cov(env0))

}
