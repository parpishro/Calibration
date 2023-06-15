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
setup_cache <- function(sim, field, Nmcmc,
                        kappa,   k0,  k1,  k2,
                        thetaS,  ts0, ts1, ts2,
                        alphaS,  as0, as1, as2,
                        thetaB,  tb0, tb1, tb2,
                        alphaB,  ab0, ab1, ab2,
                        sigma2S, ss0, ss1, ss2,
                        sigma2B, sb0, sb1, sb2,
                        sigma2E, se0, se1, se2,
                        mu,      m0,  m1,  m2) {

  cache$m   <- m <- nrow(sim)               # number of simulation runs
  cache$n   <- n <- nrow(field)             # number of field observations
  cache$p   <- p <- ncol(field) - 1         # number of experimental inputs
  cache$q   <- q <- ncol(sim) - p - 1       # number of calibration inputs
  cache$l   <- l <- q + (p + q) + (p + q) +  p + p + 1 + 1 + 1 +1 # total # of parameters


  # indices for parameters in phi
  cache$ikappa     <- ikappa   <- 1:q                                             # calibration
  cache$ithetaS    <- ithetaS  <- (q+1): (q + (p+q))                              # sim scale
  cache$ialphaS    <- ialphaS  <- (q + (p+q) + 1): (q + (p+q) + (p+q))            # sim smooth
  cache$ithetaB    <- ithetaB  <- (q + (p+q) + (p+q) + 1): (q + (p+q) + (p+q) + p)# bias scale
  cache$ialphaB    <- ialphaB  <- (q + (p+q) + (p+q) + p + 1): (l-4)              # bias smooth
  cache$isigma2S   <- isigma2S <- l-3                                             # sim variance
  cache$isigma2B   <- isigma2B <- l-2                                             # bias variance
  cache$isigma2E   <- isigma2E <- l-1                                               # error variance
  cache$imu        <- imu      <- l                                               # error variance

  # parameters (initialize first row of Phi matrix)
  Phi       <- matrix(nrow = Nmcmc, ncol = l)
  phi1      <- c(rep(k0,  length(ikappa)),
                 rep(ts0, length(ithetaS)),
                 rep(as0, length(ialphaS)),
                 rep(tb0, length(ithetaB)),
                 rep(ab0, length(ialphaB)),
                 ss0, sb0, se0, m0)
  Phi[1, ] <- phi1




  # data matrices and vectors and scaling
  Xs          <- sim[, 2:(p+q+1)]
  ys          <- sim[, 1]
  Xf          <- field[, 2:(p+1), drop=F]
  yf          <- field[, 1]

  iexp        <- 1:p
  ical        <- (p+1):(p+q)
  meanYs      <- mean(ys)
  sdYs        <- sd(ys)
  expMin      <- apply(Xs[, iexp, drop=F], 2, min)
  expMax      <- apply(Xs[, iexp, drop=F], 2, max)
  calMin      <- apply(Xs[, ical, drop=F], 2, min)
  calMax      <- apply(Xs[, ical, drop=F], 2, max)
  expRange    <- expMax-expMin
  calRange    <- calMax-calMin
  cache$scale <- list(meanYs=meanYs, sdYs=sdYs,
                      expMin=expMin, expRange=expRange,
                      calMin=calMin, calRange=calRange)

  Xs[, iexp]  <- matrix(scale(Xs[,iexp], center=expMin, scale=expRange), ncol=length(iexp))
  Xs[, ical]  <- matrix(scale(Xs[,ical], center=calMin, scale=calRange), ncol=length(ical))
  Xf          <- matrix(scale(Xf,        center=expMin, scale=expRange), ncol=length(iexp))
  y           <- as.vector(scale(c(yf, ys), center=meanYs, scale=sdYs))
  Xk          <- cbind(Xf, matrix(replicate(n, phi1[ikappa]), nrow=n, byrow=T))

  cache$Xs    <- Xs
  cache$Xk    <- Xk
  cache$Xf    <- Xf
  cache$y     <- y


  # Priors
  cache$kappa_pr   <- kappa_pr   <- setup_prior(kappa,   k1,  k2)
  cache$thetaS_pr  <- thetaS_pr  <- setup_prior(thetaS,  ts1, ts2)
  cache$alphaS_pr  <- alphaS_pr  <- setup_prior(alphaS,  as1, as2)
  cache$thetaB_pr  <- thetaB_pr  <- setup_prior(thetaB,  tb1, tb2)
  cache$alphaB_pr  <- alphaB_pr  <- setup_prior(alphaB,  ab1, ab2)
  cache$sigma2S_pr <- sigma2S_pr <- setup_prior(sigma2S, ss1, ss2)
  cache$sigma2B_pr <- sigma2B_pr <- setup_prior(sigma2B, sb1, sb2)
  cache$sigma2E_pr <- sigma2E_pr <- setup_prior(sigma2E, se1, se2)
  cache$mu_pr      <- mu_pr      <- setup_prior(mu,      m1,  m2)
  cache$priors     <- list(kappa   = c(kappa,   k1,  k2),
                           thetaS  = c(thetaS,  ts1, ts2),
                           alphaS  = c(alphaS,  as1, as2),
                           thetaB  = c(thetaB,  tb1, tb2),
                           alphaS  = c(alphaB,  ab1, ab2),
                           sigma2S = c(sigma2S, ss1, ss2),
                           sigma2B = c(sigma2B, sb1, sb2),
                           sigma2E = c(sigma2E, se1, se2),
                           mu      = c(mu,      m1,  m2))




  # Exclude parameters from MCMC
  ifixed <- c()
  if (thetaS  == "fixed") {
    Phi[, ithetaS] <- matrix(replicate(Nmcmc, phi1[ithetaS]), nrow=Nmcmc, byrow=T)
    ifixed         <- c(ifixed, ithetaS)
  }

  if (alphaS  == "fixed") {
    Phi[, ialphaS] <- matrix(replicate(Nmcmc, phi1[ialphaS]), nrow=Nmcmc, byrow=T)
    ifixed         <- c(ifixed, ialphaS)
  }

  if (thetaB  == "fixed") {
    Phi[, ithetaB] <- matrix(replicate(Nmcmc, phi1[ithetaB]), nrow=Nmcmc, byrow=T)
    ifixed         <- c(ifixed, ithetaB)
  }

  if (alphaB  == "fixed") {
    Phi[, ialphaB] <- matrix(replicate(Nmcmc, phi1[ialphaB]), nrow=Nmcmc, byrow=T)
    ifixed         <- c(ifixed, ialphaB)
  }

  if (sigma2S == "fixed")  {
    Phi[, isigma2S] <- matrix(replicate(Nmcmc, phi1[isigma2S]), nrow=Nmcmc, byrow=T)
    ifixed          <- c(ifixed, isigma2S)
  }

  if (sigma2B == "fixed")  {
    Phi[, isigma2B] <- matrix(replicate(Nmcmc, phi1[isigma2B]), nrow=Nmcmc, byrow=T)
    ifixed          <- c(ifixed, isigma2B)
  }

  if (sigma2E == "fixed")  {
    Phi[, isigma2E] <- matrix(replicate(Nmcmc, phi1[isigma2E]), nrow=Nmcmc, byrow=T)
    ifixed          <- c(ifixed, isigma2E)
  }

  if (mu == "fixed")  {
    Phi[, imu] <- matrix(replicate(Nmcmc, phi1[mu]), nrow=Nmcmc, byrow=T)
    ifixed     <- c(ifixed, imu)
  }
  cache$ifixed <- ifixed


  # set up first correlation matrices and load them into cache
  # CorKK : (n * n) correlation matrix of augmented Xk's
  # CorKS : (n * m) correlation matrix between Xk's, Xs's
  # CorSK : (m * n) correlation matrix between Xs's, Xk's
  # corSS : (m * m) correlation matrix between Xs's
  # CorFF : (n * n) correlation matrix between Xf's
  cache$CorKS  <- CorKS   <- correlation(Xk, Xs, theta = phi1[ithetaS], alpha = phi1[ialphaS])
  cache$CorKK  <- CorKK   <- correlation(Xk, Xk, theta = phi1[ithetaS], alpha = phi1[ialphaS])
  cache$CorSS  <- CorSS   <- correlation(Xs, Xs, theta = phi1[ithetaS], alpha = phi1[ialphaS])
  cache$CorFF  <- CorFF   <- correlation(Xf, Xf, theta = phi1[ithetaB], alpha = phi1[ialphaB])
  cache$CorSK  <- CorSK   <- t(CorKS)


  # compute the first log likelihood by forming the augmented covariance matrix
  #   and load both augmented covariance matirx and log likelihood into cache
  cache$Inn <- Inn <- diag(n)
  AugCov    <- rbind(cbind(phi1[isigma2S]*CorKK+phi1[isigma2B]*CorFF+phi1[isigma2E]*Inn,
                           phi1[isigma2S]*CorKS),
                     cbind(phi1[isigma2S]*CorSK, phi1[isigma2S]*CorSS))
  CholCov   <- chol(AugCov)
  InvCov    <- chol2inv(CholCov)
  logDetCov <- 2*sum(log(diag(CholCov)))




  # computes posterior log likelihood of augmented response given the augmented
  #   covariance matrix (its inverse and determinant) and residuals
  InvCovRes      <- matrix(nrow = Nmcmc, ncol = n+m)
  logPost        <- double(Nmcmc)
  res            <- y-phi1[imu]
  InvCovRes[1, ] <- InvCov%*%res
  logPost[1]     <- sum(sapply(phi1[ikappa],  kappa_pr),
                        sapply(phi1[ithetaS], thetaS_pr),
                        sapply(phi1[ialphaS], alphaS_pr),
                        sapply(phi1[ithetaB], thetaB_pr),
                        sapply(phi1[ialphaB], alphaB_pr),
                        sigma2S_pr(phi1[isigma2S]),
                        sigma2B_pr(phi1[isigma2B]),
                        sigma2E_pr(phi1[isigma2E]),
                        mu_pr(phi1[imu])) -
                    (0.5*logDetCov) -
                    (0.5*drop(t(res)%*%InvCovRes[1,]))
  cat("initial values: ", round(phi1, 3), "\n")
  return(list(Phi = Phi, logPost = logPost, InvCovRes = InvCovRes))

}




