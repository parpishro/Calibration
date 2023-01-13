update_cov <- function(covD, phi, changed, env) {

  env0             <- environment()
  parent.env(env0) <- env

  # parameter indices in phi vector



  # corFF : (n * n) correlation matrix of augmented Xf's
  # corFS : (n * m) correlation matrix between Xf's, Xs's
  # corSF : (m * n) correlation matrix between Xs's, Xf's
  # corSS : (m * m) correlation matrix between Xs's
  # corB  : (n * n) correlation matrix between Xb's

  if (changed == 0) {


  } else if (changed %in% calib) {
    covD$Xf     <- cbind(Xb, matrix(replicate(phi[calib], n), nrow = n))
    covD$CorFF  <- correlation(covD$Xf, scale = phi[scaleS], smooth = phi[smoothS])
    covD$CorFS  <- correlation(covD$Xf, Xs, scale = phi[scaleS], smooth = phi[smoothS])
    covD$CorSF  <- t(covD$CorFS)

  } else if ((changed %in% scaleS) || (changed %in% smoothS)) {
    covD$CorFF  <- correlation(covD$Xf, scale = phi[scaleS], smooth = phi[smoothS])
    covD$CorFS  <- correlation(covD$Xf, Xs, scale = phi[scaleS], smooth = phi[smoothS])
    covD$CorSF  <- t(covD$CorFS)
    covD$CorSS  <- correlation(Xs, Xs, scale = phi[scaleS], smooth = phi[smoothS])
    covD$muHat  <- mu_hat(covD$corSS, ys)
    covD$res    <- y - covD$muHat

  } else if ((changed %in% scaleB) || (changed %in% smoothB)) {
    covD$CorB   <- correlation(Xb, Xb, scale = phi[scaleB], smooth = phi[smoothB])

  } else if (changed %in% sig2S) {
    covD$sig2S  <- phi[sig2S]

  } else if (changed %in% sig2B) {
    covD$sig2B  <- phi[sig2B]

  } else if (changed %in% sigma2E) {
    covD$sig2E  <- phi[sig2E]

  } else stop("invalid changed argument!")

  In       <- diag(n)
  AugCov   <- cbind(rbind((sig2S * CorFF) + (sig2B * Xb) + (sig2E * In), CorFS),
                rbind((sig2S * CorSF), (sig2S * CorSS)))

  return(chol_cov(AugCov))

}
