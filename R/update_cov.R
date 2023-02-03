update_cov <- function(phi, changed) {
  # parameter indices in phi vector

  # indices for parameters in phi
  calib     <- 1:q
  scaleS    <- (q+1): (q + (p+q))
  smoothS   <- (q + (p+q) + 1): (q + (p+q) + (p+q))
  scaleB    <- (q + (p+q) + (p+q) + 1): (q + (p+q) + (p+q) + p)
  smoothB   <- (q + (p+q) + (p+q) + p + 1): (k-4)
  sig2S     <- k-3
  sig2B     <- k-2
  sig2E     <- k-1
  muHat     <- k #

  Xf     <- get('Xf',    envir = cache)
  Xb     <- get('Xb',    envir = cache)
  Xs     <- get('Xs',    envir = cache)
  yf     <- get('yf',    envir = cache)
  ys     <- get('ys',    envir = cache)
  y      <- get('y',     envir = cache)
  CorFF  <- get('CorFF', envir = cache)
  CorFS  <- get('CorFS', envir = cache)
  CorSF  <- get('CorSF', envir = cache)
  CorSS  <- get('CorSS', envir = cache)
  CorB   <- get('CorB',  envir = cache)
  muHat  <- get('muHat', envir = cache)
  res    <- get('res',   envir = cache)
  res    <- get('res',   envir = cache)






  # corFF : (n * n) correlation matrix of augmented Xf's
  # corFS : (n * m) correlation matrix between Xf's, Xs's
  # corSF : (m * n) correlation matrix between Xs's, Xf's
  # corSS : (m * m) correlation matrix between Xs's
  # corB  : (n * n) correlation matrix between Xb's

  if (changed %in% calib) {
    Xf     <- cbind(Xb, replicate(n, phi[calib], n))
    CorFF  <- corelation(Xf, scale = phi[scaleS], smooth = phi[smoothS])
    CorFS  <- corelation(Xf, Xs, scale = phi[scaleS], smooth = phi[smoothS])
    CorSF  <- t(CorFS)

    assign('Xf', Xf, envir = cache)
    assign('CorFF', CorFF, envir = cache)
    assign('CorFS', CorFS, envir = cache)
    assign('CorSF', CorSF, envir = cache)

  } else if ((changed %in% scaleS) || (changed %in% smoothS)) {
    CorFF  <- corelation(Xf, scale = phi[scaleS], smooth = phi[smoothS])
    CorFS  <- corelation(Xf, Xs, scale = phi[scaleS], smooth = phi[smoothS])
    CorSF  <- t(CorFS)
    CorSS  <- corelation(Xs, scale = phi[scaleS], smooth = phi[smoothS])
    muHat  <- mu_hat(CorSS, ys)
    res    <- y - muHat

    assign('CorSS', CorSS, envir = cache)
    assign('CorFF', CorFF, envir = cache)
    assign('CorFS', CorFS, envir = cache)
    assign('CorSF', CorSF, envir = cache)
    assign('muHat', muHat, envir = cache)
    assign('res',   res,   envir = cache)

  } else if ((changed %in% scaleB) || (changed %in% smoothB)) {
    CorB   <- corelation(Xb, scale = phi[scaleB], smooth = phi[smoothB])
    assign('CorB', CorB, envir = cache)

  } else if (changed %in% sig2S) {
    sig2S  <- phi[sig2S]
    assign('sig2S', sig2S, envir = cache)

  } else if (changed %in% sig2B) {
    sig2B  <- phi[sig2B]
    assign('sig2B', sig2B, envir = cache)

  } else if (changed %in% sig2E) {
    sig2E  <- phi[sig2E]
    assign('sig2E', sig2E, envir = cache)


  } else stop("invalid changed argument!")

  I        <- diag(n)
  AugCov   <- rbind(cbind(((sig2S * CorFF) + (sig2B * CorB) + (sig2E * I)), (sig2S * CorFS)),
                    cbind((sig2S * CorSF), (sig2S * CorSS)))
  CholCov   <- chol(AugCov)
  InvCov    <- chol2inv(CholCov)
  logDetCov <- sum(2*log(diag(CholCov)))

  return(list(InvCov = InvCov, logDetCov = logDetCov))

}
