#' Sets up the cache environment to hold all the intermediate results
#'
#' `setup_cache` declares and initalize the basic information about the data such as
#' number of parameters in each category. It sets up all the initial values of the parameters
#' and runs a first iteration of the workflow to compute log posterior and all the intermediate
#' objects such as correlation matrices of both GPs. The objects in cache environment can be
#' modified/accessed from all package functions.
#'
#' @param sim       a matrix containing the simulation data
#' @param field     a matrix containing the field data
#' @param kappaPr   a function that computes log of prior for calibration parameters
#' @param thetaPr   a function that computes log of prior for scale hyperparameters
#' @param alphaPr   a function that computes log of prior for smoothness hyperparameters
#' @param sigma2Pr  a function that computes log of prior for variance hyperparameters
#' @noRd
#' @return A list consisting of sampled Parameters in matrix form and log posterior vector
setup_cache <- function(sim, field, kappaPr,thetaPr, alphaPr, sigma2Pr) {

  m   <- nrow(sim)               # number of simulation runs
  n   <- nrow(field)             # number of field observations
  p   <- ncol(field) - 1         # number of experimental variables
  q   <- ncol(sim) - p - 1       # number of calibration parameters
  d   <- ncol(sim) - 1           # number of all variables for simulation
  k   <- q + (p + q) + (p + q) +  p + p + 1 + 1 + 1 # total # of parameters

  assign('m', m, envir = cache)    # number of simulation runs
  assign('n', n, envir = cache)    # number of field observations
  assign('p', p, envir = cache)    # number of experimental variables
  assign('q', q, envir = cache)    # number of calibration parameters
  assign('d', d, envir = cache)    # number of all variables for simulation
  assign('k', k, envir = cache)

  # indices for parameters in phi
  ikappa     <- 1:q                                             # calibration
  ithetaS    <- (q+1): (q + (p+q))                              # sim scale
  ialphaS    <- (q + (p+q) + 1): (q + (p+q) + (p+q))            # sim smooth
  ithetaB    <- (q + (p+q) + (p+q) + 1): (q + (p+q) + (p+q) + p)# bias scale
  ialphaB    <- (q + (p+q) + (p+q) + p + 1): (k-3)              # bias smooth
  isigma2S   <- k-2                                             # sim variance
  isigma2B   <- k-1                                             # bias variance
  isigma2E   <- k                                             # error variance


  assign('ikappa',   ikappa ,  envir = cache)     # calibration
  assign('ithetaS',  ithetaS,  envir = cache)     # sim scale
  assign('ialphaS',  ialphaS,  envir = cache)     # sim smoothness
  assign('ithetaB',  ithetaB,  envir = cache)     # bias scale
  assign('ialphaB',  ialphaB,  envir = cache)     # bias smoothness
  assign('isigma2S', isigma2S, envir = cache)     # sim variance
  assign('isigma2B', isigma2B, envir = cache)     # bias variance
  assign('isigma2E', isigma2E, envir = cache)     # error variance
  #assign('imuHat',   imuHat,   envir = cache)     # number of total parameters

  # data matrices and vectors and scaling
  Xs      <- sim[, 1:d]
  ys      <- sim[, d + 1]
  ys      <- sim[, d + 1]
  Xb      <- field[, 1:p, drop = FALSE]
  yf      <- field[, p + 1]

  y       <- scale(matrix(c(yf, ys), ncol = 1), center =  mean(ys), scale = sd(ys))

  calibMin      <- apply(Xs[, (p+1):d, drop=F], 2, min)
  calibMax      <- apply(Xs[, (p+1):d, drop=F], 2, max)
  Xs[, (p+1):d] <- scale(Xs[, (p+1):d, drop=F], center = calibMin, scale = calibMax - calibMin)

  #             <- rbind(Xs[, 1:p, drop=F], Xb)
  experMin      <- apply(Xs[, 1:p, drop=F], 2, min)
  experMax      <- apply(Xs[, 1:p, drop=F], 2, max)
  Xs[, 1:p]     <- scale(Xs[, 1:p, drop=F], center = experMin, scale = experMax - experMin)
  Xb            <- scale(Xb,                center = experMin, scale = experMax - experMin)

  assign('calibMin', calibMin, envir = cache)
  assign('calibMax', calibMax, envir = cache)
  assign('Xs',       Xs,       envir = cache)
  assign('ys',       ys,       envir = cache)
  assign('Xb',       Xb,       envir = cache)
  assign('yf',       yf,       envir = cache)
  assign('y',        y,        envir = cache)

  # parameters (initialize first row of Phi matrix)
  phi            <- double(k)
  phi[ikappa]    <- apply(Xs[, (p+1):(p+q), drop=F], 2, mean)  #0.5
  phi[ithetaS]   <- double(length(ithetaS)) + 0.5
  phi[ialphaS]   <- double(length(ialphaS)) + 1.6
  phi[ithetaB]   <- double(length(ithetaB)) + 4
  phi[ialphaB]   <- double(length(ialphaB)) + 1.86
  phi[isigma2S]  <- 3.6
  phi[isigma2B]  <- 1.6
  phi[isigma2E]  <- 0.1

  # set up first correlation matrices and load them into cache
  # corFF : (n * n) correlation matrix of augmented Xf's
  # corFS : (n * m) correlation matrix between Xf's, Xs's
  # corSF : (m * n) correlation matrix between Xs's, Xf's
  # corSS : (m * m) correlation matrix between Xs's
  # corB  : (n * n) correlation matrix between Xb's
  Xf      <- cbind(Xb, replicate(n, phi[ikappa], n))
  CorFS   <- correlation(Xf, Xs, theta = phi[ithetaS], alpha = phi[ialphaS])
  CorFF   <- correlation(Xf,     theta = phi[ithetaS], alpha = phi[ialphaS])
  CorSS   <- correlation(Xs,     theta = phi[ithetaS], alpha = phi[ialphaS])
  CorB    <- correlation(Xb,     theta = phi[ithetaB], alpha = phi[ialphaB])
  CorSF   <- t(CorFS)
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
  AugCov      <- rbind(cbind(sigma2S*CorFF + sigma2B*CorB + sigma2E*Inn, sigma2S*CorFS),
                       cbind((sigma2S*CorSF), (sigma2S*CorSS)))
  CholCov     <- chol(AugCov)
  InvCov      <- chol2inv(CholCov)
  logDetCov   <- 2*sum(log(diag(CholCov)))

  assign('InvCov',  InvCov,  envir = cache)
  #muHat       <- update_mu()
  #phi[imuHat] <- muHat

  #assign('muHat',   muHat,   envir = cache)

  # computes posterior log likelihood of augmented response given the augmented
  #   covariance matrix (its inverse and determinant) and residuals
  lPost <- sum(sapply(phi[ikappa],              kappaPr)  +
               sapply(phi[c(ithetaS, ithetaB)], thetaPr)  +
               sapply(phi[c(ialphaS, ialphaB)], alphaPr)  +
               sapply(phi[isigma2S:isigma2E],   sigma2Pr)) - (0.5*logDetCov + 0.5 * drop(t(y)%*%InvCov%*%y))
  return(list(phi = phi, logPost = lPost))

}




