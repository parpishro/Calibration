#' setup_cache
#' sets up all the initial values based on given data and priors
#'
#' @param sim
#' @param field
#' @param thetaPr
#' @param lambdaPr
#' @param gammaPr
#' @param sigma2Pr
#'
#' @return
setup_cache <- function(sim, field, thetaPr,lambdaPr, gammaPr, sigma2Pr) {

  # setting cache environment: TO BE ACCESSED/MODIFIED BY ALL PACKAGE FUNCTIONS
  # scalers
  m   <- nrow(sim)               # number of simulation runs
  n   <- nrow(field)             # number of field observations
  p   <- ncol(field) - 1         # number of experimental variables
  q   <- ncol(sim) - p - 1       # number of calibration parameters
  d   <- ncol(sim) - 1           # number of all variables for simulation
  k   <- q + (p + q) + (p + q) +  p + p + 1 + 1 + 1 + 1 # total # of parameters


  assign('m', m, envir = cache)    # number of simulation runs
  assign('n', n, envir = cache)    # number of field observations
  assign('p', p, envir = cache)    # number of experimental variables
  assign('q', q, envir = cache)    # number of calibration parameters
  assign('d', d, envir = cache)    # number of all variables for simulation
  assign('k', k, envir = cache)

  # indices for parameters in phi
  itheta     <- 1:q                                             # calibration
  ilambdaS   <- (q+1): (q + (p+q))                              # sim scale
  igammaS    <- (q + (p+q) + 1): (q + (p+q) + (p+q))            # sim smooth
  ilambdaB   <- (q + (p+q) + (p+q) + 1): (q + (p+q) + (p+q) + p)# bias scale
  igammaB    <- (q + (p+q) + (p+q) + p + 1): (k-4)              # bias smooth
  isigma2S   <- k-3                                             # sim variance
  isigma2B   <- k-2                                             # bias variance
  isigma2E   <- k-1                                             # error variance
  imuHat     <- k                                               # estimated mean

  assign('itheta',   itheta ,  envir = cache)     # calibration
  assign('ilambdaS', ilambdaS, envir = cache)     # sim scale
  assign('igammaS',  igammaS,  envir = cache)     # sim smoothness
  assign('ilambdaB', ilambdaB, envir = cache)     # bias scale
  assign('igammaB',  igammaB,  envir = cache)     # bias smoothness
  assign('isigma2S', isigma2S, envir = cache)     # sim variance
  assign('isigma2B', isigma2B, envir = cache)     # bias variance
  assign('isigma2E', isigma2E, envir = cache)     # error variance
  assign('imuHat',   imuHat,   envir = cache)     # number of total parameters

  # data matrices and vectors
  Xs      <- sim[, 1:d]
  ys      <- sim[, d + 1]
  Xb      <- field[, 1:p, drop = FALSE]
  yf      <- field[, p + 1]
  y       <- (c(ys, yf) - mean(ys)) / sd(ys)

  assign('Xs', Xs, envir = cache)
  assign('ys', ys, envir = cache)
  assign('Xb', Xb, envir = cache)
  assign('yf', yf, envir = cache)
  assign('y',  y,  envir = cache)


  # parameters (initialize first row of Phi matrix)
  phi            <- double(k)
  phi[itheta]    <- double(length(itheta))  + thetaPr$mean
  phi[ilambdaS]  <- double(length(ilambdaS)) + lambdaPr$mean
  phi[igammaS]   <- double(length(igammaS)) + gammaPr$mean
  phi[ilambdaB]  <- double(length(ilambdaB)) + lambdaPr$mean
  phi[igammaB]   <- double(length(igammaB)) + gammaPr$mean
  phi[isigma2S]  <- sigma2Pr$mean
  phi[isigma2B]  <- sigma2Pr$mean
  phi[isigma2E]  <- sigma2Pr$mean



  # set up first correlation matrices and load them into cache
  Xf      <- cbind(Xb, replicate(n, phi[itheta], n))
  CorFF   <- correlation(Xf,     lambda = phi[ilambdaS], gamma = phi[igammaS])
  CorFS   <- correlation(Xf, Xs, lambda = phi[ilambdaS], gamma = phi[igammaS])
  CorSF   <- t(CorFS)
  CorSS   <- correlation(Xs,     lambda = phi[ilambdaS], gamma = phi[igammaS])
  CorB    <- correlation(Xb,     lambda = phi[ilambdaB], gamma = phi[igammaB])


  sigma2S <- phi[isigma2S]
  sigma2B <- phi[isigma2B]
  sigma2E <- phi[isigma2E]


  assign('Xf',      Xf,      envir = cache)
  assign('CorFF',   CorFF,   envir = cache)
  assign('CorFS',   CorFS,   envir = cache)
  assign('CorSF',   CorSF,   envir = cache)
  assign('CorSS',   CorSS,   envir = cache)
  assign('CorB',    CorB,    envir = cache)
  assign('sigma2S', sigma2S, envir = cache)
  assign('sigma2B', sigma2B, envir = cache)
  assign('sigma2E', sigma2E, envir = cache)


  # compute the first log likelihood by forming the augmented covariance matrix
  #   and load both augmented covariance matirx and log likelihood into cache
  Inn         <- diag(n)
  AugCov      <- rbind(cbind(sigma2S*CorFF + sigma2B*CorB + sigma2E*Inn,
                             sigma2S*CorFS),
                       cbind((sigma2S*CorSF), (sigma2S*CorSS)))
  CholCov     <- chol(AugCov)
  InvCov      <- chol2inv(CholCov)
  logDetCov   <- sum(2*log(diag(CholCov)))

  assign('Chol',  CholCov,  envir = cache)
  muHat       <- update_mu()
  res         <- y - muHat
  phi[imuHat] <- muHat

  assign('muHat',   muHat,   envir = cache)
  assign('res',     res,     envir = cache)



  # computes posterior log likelihood of augmented response given the augmented
  #   covariance matrix (its inverse and determinant) and residuals
  lPost <- sum(log(phi)) - (0.5*logDetCov) - (0.5 * (res%*%InvCov%*%res))

  return(list(phi = phi, logPost = lPost))

}
