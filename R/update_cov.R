update_cov <- function(covD, phi, changed, env) {

  env0             <- environment()
  parent.env(env0) <- env

  # parameter indices in phi vector
  calib    <- 1:q
  scaleS   <- (q+1): (q + (p + q))
  smoothS  <- (q + (p + q) + 1): (q + (p + q) + (p + q))
  scaleB   <- (q + (p + q) + (p + q) + 1): (q + (p + q) + (p + q) + p)
  smoothB  <- (q + (p + q) + (p + q) + p + 1): (k - 4)
  sig2S    <- k - 3
  sig2B    <- k - 2
  sig2E    <- k - 1
  muHat    <- k


  # corFF : (n * n) correlation matrix of augmented Xf's
  # corFS : (n * m) correlation matrix between Xf's, Xs's
  # corSF : (m * n) correlation matrix between Xs's, Xf's
  # corSS : (m * m) correlation matrix between Xs's
  # corB  : (n * n) correlation matrix between Xb's
  Xf       <- covD$Xf
  CorFF    <- covD$CorFF
  CorFS    <- covD$CorFS
  CorSF    <- covD$CorSF
  CorSS    <- covD$CorSS
  muHat    <- covD$muHat
  res      <- covD$res
  CorB     <- covD$CorB

  if (changed == 0) {
    Xf     <- cbind(Xb, matrix(replicate(phi[calib], n), nrow = n))
    CorFF  <- correlation(Xf, scale = phi[scaleS], smooth = phi[smoothS])
    CorFS  <- correlation(Xf, Xs, scale = phi[scaleS], smooth = phi[smoothS])
    CorSF  <- t(CorFS)
    CorSS  <- correlation(Xs, Xs, scale = phi[scaleS], smooth = phi[smoothS])
    muHat  <- mu_hat(corSS, ys)
    res    <- z - muHat
    CorB   <- correlation(Xb, Xb, scale = phi[scaleB], smooth = phi[smoothB])

  } else if (changed %in% calib) {
    Xf     <- cbind(Xb, matrix(replicate(phi[calib], n), nrow = n))
    CorFF  <- correlation(Xf, Xf, scale = phi[scaleS], smooth = phi[smoothS])
    CorFS  <- correlation(Xf, Xs, scale = phi[scaleS], smooth = phi[smoothS])
    CorSF  <- t(CorFS)

  } else if ((changed %in% scaleS) || (changed %in% smoothS)) {
    CorFF  <- correlation(Xf, Xf, scale = phi[scaleS], smooth = phi[smoothS])
    CorFS  <- correlation(Xf, Xs, scale = phi[scaleS], smooth = phi[smoothS])
    CorSF  <- t(CorFS)
    CorSS  <- correlation(Xs, Xs, scale = phi[scaleS], smooth = phi[smoothS])
    muHat  <- mu_hat(corSS, ys)
    res    <- z - muHat

  } else if ((changed %in% scaleB) || (changed %in% smoothB)) {
    CorB   <- correlation(Xb, Xb, scale = phi[scaleB], smooth = phi[smoothB])

  } else if (changed %in% sig2S) {
    sig2S  <- phi[sig2S]

  } else if (changed %in% sig2B) {
    sig2B  <- phi[sig2B]

  } else if (changed %in% sigma2E) {
    sig2E  <- phi[sig2E]

  } else stop("invalid changed argument!")

  In.      <- diag(n)
  AugCov   <- cbind(rbind((sig2S * CorFF) + (sig2B * Xb) + (sig2E * In), CorFS),
                rbind((sig2S * CorSF), (sig2S * CorSS)))
  chol     <- chol_cov(AugCov)
  inv <- chol$inv
  det <- chol$det

  return(list(Xf = Xf, CorFF = CorFF, CorFS = CorFS, CorSF = CorSf,
              CorSS = CorSS, muHat = muHat, res = res,
              CorB = CorB, sig2S = sig2S, sig2B = sig2B, sig2E = sig2E))

}
