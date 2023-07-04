#' Setting Cache Environment
#'
#' `setup_cache()` sets up a cache environment to hold shared data. The cache can be
#' accessed and modified by all package function. Any large dataset appearing iteratively
#' must be kept here to minimize argument passing. Furthermore, it sets up the prior
#' functions and initiates the Phi matrix to be completed by MCMC algorithm.
#'
#' @param sim       numeric matrix containing the simulation data
#' @param field     numeric matrix containing the field data
#' @param priors    nested list containing prior specification for all parameters
#' @param Nmcmc     integer representing number of MCMC runs
#'
#' @return A list consisting of initialized `Phi` matrix and `logPost` vector, and prior
#' function list `priorFns`
#' @noRd
setup_cache <- function(sim, field, priors, Nmcmc) {

  cache$m   <- m <- nrow(sim)               # number of simulation runs
  cache$n   <- n <- nrow(field)             # number of field observations
  cache$p   <- p <- ncol(field) - 1         # number of experimental inputs
  cache$q   <- q <- ncol(sim) - p - 1       # number of calibration inputs
  cache$l   <- l <- q + (p+q) + (p+q) +  p + p + 1+1+1+1 # total # of parameters

  # indices for parameters in phi
  cache$ikappa   <- ikappa   <- 1:q                                             # calibration
  cache$ithetaS  <- ithetaS  <- (q+1): (q + (p+q))                              # sim scale
  cache$ialphaS  <- ialphaS  <- (q + (p+q) + 1): (q + (p+q) + (p+q))            # sim smooth
  cache$ithetaB  <- ithetaB  <- (q + (p+q) + (p+q) + 1): (q + (p+q) + (p+q) + p)# bias scale
  cache$ialphaB  <- ialphaB  <- (q + (p+q) + (p+q) + p + 1): (l-4)              # bias smooth
  cache$isigma2S <- isigma2S <- l-3                                             # sim variance
  cache$isigma2B <- isigma2B <- l-2                                             # bias variance
  cache$isigma2E <- isigma2E <- l-1                                             # error variance
  cache$imuB     <- imuB     <- l                                               # error variance
  indices        <- list(ikappa   =ikappa,
                         ithetaS  = ithetaS, ialphaS = ialphaS,
                         ithetaB  = ithetaB, ialphaB = ialphaB,
                         isigma2S = isigma2S, isigma2B = isigma2B,
                         isigma2E = isigma2E, imuB     = imuB)
  cache$indices  <- indices



  # parameters (initialize first row of Phi matrix)
  Phi      <- matrix(nrow = Nmcmc, ncol = l)
  ifixed   <- c()
  phi1     <- c()
  priorFns <- list()
  for (i in 1:length(priors)) {
    param    <- priors[[i]]
    iparam   <- indices[[i]]
    dist <- param$dist
    init <- param$init
    p1   <- param$p1
    p2   <- param$p2

    if (length(init) == 1)
      phi1     <- c(phi1, rep(init, length(iparam)))
    else if (length(init) == length(iparam))
      phi1     <- c(phi1, init)
    else
      stop("Length of parameter initial values must be same as parameter elements!")

    if (dist  == "fixed") {
      Phi[, iparam] <- matrix(rep(phi1[iparam], Nmcmc), nrow=Nmcmc, byrow=T)
      ifixed        <- c(ifixed, iparam)
      priorFns[iparam] <-  replicate(length(iparam), function(x) 0)
    } else if (length(dist) ==  1) {
      priorFns[iparam] <- replicate(length(iparam), prior_builder(dist, p1, p2))
    } else if (length(dist) ==  length(iparam)) {
      for (j in 1:length(dist)) {
        priorFns    <- c(priorFns, prior_builder(dist[j], p1, p2))
      }

    } else
        stop("Length of parameter priors distributions must be same as parameter elements!")
  }

  cache$ifixed   <- ifixed
  cache$priorFns <- priorFns
  cache$priors   <- priors


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
  cache$scale <- list(meanYs=meanYs, sdYs=sdYs, expMin=expMin, expRange=expRange, calMin=calMin, calRange=calRange)
  phi1[ikappa]<- (phi1[ikappa] - calMin) / calRange
  Xs[, iexp]  <- matrix(scale(Xs[,iexp], center=expMin, scale=expRange), ncol=length(iexp))
  Xs[, ical]  <- matrix(scale(Xs[,ical], center=calMin, scale=calRange), ncol=length(ical))
  Xf          <- matrix(scale(Xf,        center=expMin, scale=expRange), ncol=length(iexp))
  y           <- as.vector(scale(c(yf, ys), center=meanYs, scale=sdYs))

  TrueKappa   <- matrix(replicate(n, phi1[ikappa]), nrow=n, byrow=T)
  Xk          <- cbind(Xf, TrueKappa)

  cache$Xs    <- Xs
  cache$Xk    <- Xk
  cache$Xf    <- Xf
  cache$y     <- y

  # proposal sd rates
  cache$sdRates <- c(rep(0.2, length(ikappa)), rep(0.2, length(ithetaS)), rep(0.2, length(ialphaS)),
                     rep(0.2, length(ithetaB)), rep(0.2, length(ialphaB)), 0.2, 0.2, 0.2, 0.2)

  # set up first correlation matrices and load them into cache
  # CorKK : (n * n) correlation matrix of augmented Xk's
  # CorKS : (n * m) correlation matrix between Xk's, Xs's
  # CorSK : (m * n) correlation matrix between Xs's, Xk's
  # corSS : (m * m) correlation matrix between Xs's
  # CorFF : (n * n) correlation matrix between Xf's
  cache$CorKK  <- CorKK   <- correlation(Xk, Xk, theta = phi1[ithetaS], alpha = phi1[ialphaS])
  cache$CorKS  <- CorKS   <- correlation(Xk, Xs, theta = phi1[ithetaS], alpha = phi1[ialphaS])
  cache$CorSK  <- CorSK   <- t(CorKS)
  cache$CorSS  <- CorSS   <- correlation(Xs, Xs, theta = phi1[ithetaS], alpha = phi1[ialphaS])
  cache$CorFF  <- CorFF   <- correlation(Xf, Xf, theta = phi1[ithetaB], alpha = phi1[ialphaB])

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
  logPost        <- double(Nmcmc)
  res            <- y-phi1[imuB]
  logPrior       <- 0
  for (i in 1:length(priorFns)) {
    fn <- priorFns[[i]]
    logPrior     <- logPrior + fn(phi1[[i]])
  }
  logPost[1]     <- logPrior - 0.5*(logDetCov+drop(t(res)%*%InvCov%*%res))
  Phi[1, ]       <- phi1
  cat("initial values: ", round(phi1, 3), "\n")
  return(list(Phi = Phi, logPost = logPost))
}




