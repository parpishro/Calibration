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
#' @param k0        a scaler or vector of same length as number of calibration
#'                  parameters that represent initial value for calibration parameters
#'                  after scaling to (0, 1)
#' @param t0        a scaler in (0, Inf) representing initial values of theta (scale) parameters
#' @param a0        a scaler in (1, 2) representing initial values of alpha (smoothness) parameters
#' @param s0        a scaler or vector of length three that represent initial value for marginal variance
#'
#' @param sigma2Pr  a function that computes log of prior for variance hyperparameters
#'
#' @noRd
#' @return A list consisting of sampled Parameters in matrix form and log posterior vector
setup_cache <- function(sim, field, kappaPr,thetaPr, alphaPr, sigma2Pr, k0, t0, a0, s0) {

  cache$m   <- m <- nrow(sim)               # number of simulation runs
  cache$n   <- n <- nrow(field)             # number of field observations
  cache$p   <- p <- ncol(field) - 1         # number of experimental variables
  cache$q   <- q <- ncol(sim) - p - 1       # number of calibration parameters
  cache$d   <- d <- ncol(sim) - 1           # number of all variables for simulation
  cache$k   <- k <- q + (p + q) + (p + q) +  p + p + 1 + 1 + 1 # total # of parameters


  # indices for parameters in phi
  cache$ikappa     <- ikappa   <- 1:q                                             # calibration
  cache$ithetaS    <- ithetaS  <- (q+1): (q + (p+q))                              # sim scale
  cache$ialphaS    <- ialphaS  <- (q + (p+q) + 1): (q + (p+q) + (p+q))            # sim smooth
  cache$ithetaB    <- ithetaB  <- (q + (p+q) + (p+q) + 1): (q + (p+q) + (p+q) + p)# bias scale
  cache$ialphaB    <- ialphaB  <- (q + (p+q) + (p+q) + p + 1): (k-3)              # bias smooth
  cache$isigma2S   <- isigma2S <- k-2                                             # sim variance
  cache$isigma2B   <- isigma2B <- k-1                                             # bias variance
  cache$isigma2E   <- isigma2E <- k                                             # error variance


  # data matrices and vectors and scaling
  cache$Xs            <- Xs         <- sim[, 1:d]
  cache$ys            <- ys         <- sim[, d+1]
  cache$Xb            <- Xb         <- field[, 1:p, drop=F]
  cache$yf            <- yf         <- field[, p+1]

  cache$y             <- y          <- scale(matrix(c(yf, ys), ncol=1), center=mean(ys), scale=sd(ys))

  cache$calibMin      <- calibMin   <- apply(Xs[, (p+1):d, drop=F], 2, min)
  cache$calibMax      <- calibMax   <- apply(Xs[, (p+1):d, drop=F], 2, max)
  cache$Xs[, (p+1):d] <- scale(Xs[, (p+1):d, drop=F], center=calibMin, scale=calibMax-calibMin)

  cache$experMin      <- experMin   <- apply(Xs[, 1:p, drop=F], 2, min)
  cache$experMax      <- experMax   <- apply(Xs[, 1:p, drop=F], 2, max)
  cache$Xs[, 1:p]     <- scale(Xs[, 1:p, drop=F], center=experMin, scale=experMax-experMin)
  cache$Xb            <- Xb         <- scale(Xb, center=experMin, scale=experMax-experMin)
  Xs                  <- cache$Xs



  # parameters (initialize first row of Phi matrix)
  phi            <- double(k)
  phi[ikappa]    <- double(length(ikappa)) + k0
  phi[ithetaS]   <- double(length(ithetaS)) + t0
  phi[ialphaS]   <- double(length(ialphaS)) + a0
  phi[ithetaB]   <- double(length(ithetaB)) + t0
  phi[ialphaB]   <- double(length(ialphaB)) + a0
  phi[isigma2S:isigma2E]  <- double(length(isigma2S:isigma2E)) + s0




  # set up first correlation matrices and load them into cache
  # corFF : (n * n) correlation matrix of augmented Xf's
  # corFS : (n * m) correlation matrix between Xf's, Xs's
  # corSF : (m * n) correlation matrix between Xs's, Xf's
  # corSS : (m * m) correlation matrix between Xs's
  # corB  : (n * n) correlation matrix between Xb's
  cache$Xf      <- Xf      <- cbind(Xb, replicate(n, phi[ikappa], n))
  cache$CorFS   <- CorFS   <- correlation(Xf, Xs, theta = phi[ithetaS], alpha = phi[ialphaS])
  cache$CorFF   <- CorFF   <- correlation(Xf,     theta = phi[ithetaS], alpha = phi[ialphaS])
  cache$CorSS   <- CorSS   <- correlation(Xs,     theta = phi[ithetaS], alpha = phi[ialphaS])
  cache$CorB    <- CorB    <- correlation(Xb,     theta = phi[ithetaB], alpha = phi[ialphaB])
  cache$CorSF   <- CorSF   <- t(CorFS)
  cache$sigma2S <- sigma2S <- phi[isigma2S]
  cache$sigma2B <- sigma2B <- phi[isigma2B]
  cache$sigma2E <- sigma2E <- phi[isigma2E]


  # compute the first log likelihood by forming the augmented covariance matrix
  #   and load both augmented covariance matirx and log likelihood into cache
  Inn             <- diag(n)
  AugCov          <- rbind(cbind(sigma2S*CorFF + sigma2B*CorB + sigma2E*Inn, sigma2S*CorFS),
                           cbind((sigma2S*CorSF), (sigma2S*CorSS)))
  CholCov         <- chol(AugCov)
  InvCov    <- chol2inv(CholCov)
  logDetCov <- 2*sum(log(diag(CholCov)))




  # computes posterior log likelihood of augmented response given the augmented
  #   covariance matrix (its inverse and determinant) and residuals
  lPost <- sum(sapply(phi[ikappa],              kappaPr)  +
               sapply(phi[c(ithetaS, ithetaB)], thetaPr)  +
               sapply(phi[c(ialphaS, ialphaB)], alphaPr)  +
               sapply(phi[isigma2S:isigma2E],   sigma2Pr)) - (0.5*logDetCov + 0.5 * drop(t(y)%*%InvCov%*%y))
  return(list(phi = phi, logPost = lPost))

}




