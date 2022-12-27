update_cov <- function(covD, phi, changed, env) {

  env0             <- environment()
  parent.env(env0) <- env

  # parameter indices in phi vector
  calib            <- 1:q
  scaleS           <- (q+1): (q + (p + q))
  smoothS          <- (q + (p + q) + 1): (q + (p + q) + (p + q))
  scaleB           <- (q + (p + q) + (p + q) + 1): (q + (p + q) + (p + q) + p)
  smoothB          <- (q + (p + q) + (p + q) + p + 1): (k - 4)
  sigma2S          <- k - 3
  sigma2B          <- k - 2
  sigma2E          <- k - 1
  muHat            <- k




  corFF            <- covD$corFF
  corFS            <- covD$corFS
  corSF            <- covD$corSF
  corSS            <- covD$corSS
  corB             <- covD$corB


  if (changed == 0) {

    # (n * n) correlation matrix between augmented Xf's
    Xf     <- cbind(Xb, matrix(replicate(phi[calib], n), nrow = n))
    CorFF  <- correlation(Xf, scale = phi[scaleS], smooth = phi[smoothS])

    # (n * m) correlation matrix between Xf's, Xs's
    CorFS  <- correlation(Xf, Xs, , scale = phi[scaleS], smooth = phi[smoothS])

    # (m * n) correlation matrix between Xs's, Xf's
    CorSF  <- t(CorFS)

    # (m * m) correlation matrix between Xs's
    CorSS  <- correlation(Xs, scale = phi[scaleS], smooth = phi[smoothS])

    # (n * n) correlation matrix between Xb's
    CorB   <- correlation(Xb, scale = phi[scaleB], smooth = phi[smoothB])

  } else if (changed == "omega sim") {

  } else if (changed == "alpha sim") {

  } else if (changed == "omega bias") {

  } else if (changed == "alpha bias") {

  } else if (changed == "variance sim") {

  } else if (changed == "variance bias") {

  } else if (changed == "variance sim") {

  } else if (is.null(changed)) {

  } else stop("invalid parameter type!")

}
