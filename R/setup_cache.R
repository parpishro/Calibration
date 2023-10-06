#' Setting Cache Environment
#'
#' `setup_cache()` sets up a cache environment to hold shared data. The cache can be
#' accessed and modified by all package function. Any large dataset appearing iteratively
#' must be kept here to minimize argument passing. Furthermore, it sets up the prior
#' functions and initiates the Phi matrix to be completed by MCMC algorithm.
#'
#' @param sim       numeric matrix containing the simulation data
#' @param field     numeric matrix containing the field data
#' @param priorFamilies    nested list containing prior specification for all parameters
#' @param nMCMC     integer representing number of MCMC runs
#' @param showProgress  logical indicating whether progress must be displayed at console.
#'                    Default is False.
#'
#' @return A list consisting of initialized `Phi` matrix and `logPost` vector
#' @importFrom stats sd
#' @noRd
setup_cache <- function(sim, field, priorFamilies, nMCMC, showProgress) {

  cache$m   <- m <- nrow(sim)               # number of simulation runs
  cache$n   <- n <- nrow(field)             # number of field observations
  cache$p   <- p <- ncol(field) - 1         # number of experimental inputs
  cache$q   <- q <- ncol(sim) - p - 1       # number of calibration inputs
  cache$l   <- l <- q + (p + q) + (p + q) +  p + p + 1 + 1 + 1 + 1 # total # of parameters

  # indices for parameters in phi
  cache$ikappa   <- ikappa   <- 1:q                                                    # calibration
  cache$ithetaS  <- ithetaS  <- (q + 1):(q + (p + q))                                 # sim scale
  cache$ialphaS  <- ialphaS  <- (q + (p + q) + 1):(q + (p + q) + (p + q))              # sim smooth
  cache$ithetaB  <- ithetaB  <- (q + (p + q) + (p + q) + 1):(q + (p + q) + (p + q) + p)# bias scale
  cache$ialphaB  <- ialphaB  <- (q + (p + q) + (p + q) + p + 1):(l - 4)                # bias smooth
  cache$isigma2S <- isigma2S <- l - 3                                                    # sim variance
  cache$isigma2B <- isigma2B <- l - 2                                                    # bias variance
  cache$isigma2E <- isigma2E <- l - 1                                                    # error variance
  cache$imuB     <- imuB     <- l                                               # error variance
  cache$iexp     <- iexp     <- 1:p
  cache$ical     <- ical     <- (p + 1):(p + q)
  cache$indices  <- indices  <- list(ikappa   = ikappa,
                                     ithetaS  = ithetaS,  ialphaS  = ialphaS,
                                     ithetaB  = ithetaB,  ialphaB  = ialphaB,
                                     isigma2S = isigma2S, isigma2B = isigma2B,
                                     isigma2E = isigma2E, imuB     = imuB,
                                     n = n, m = m, p = p, q = q, l = l)


  # setting up priors based on user-given specification
  Phi      <- matrix(nrow = nMCMC, ncol = l)
  ifixed   <- c()
  priors   <- list(param = c(), dist = c(), p1  = c(),  p2 = c(),
                   mean  = c(), sd   = c(), fun = list())

  for (i in 1:length(priorFamilies)) {
    iparam   <- indices[[i]]
    priorFam <- priorFamilies[[i]]
    dist     <- priorFam$dist
    p1       <- priorFam$p1
    p2       <- priorFam$p2
    for (j in iparam) {
      ind             <- j - iparam[1] + 1
      priors$param[j] <- if (i < 6) paste0(names(priorFamilies)[i], ind)
                         else  names(priorFamilies)[i]
      priors$dist[j]  <- if (length(dist) == 1) dist else dist[ind]
      priors$p1[j]    <- if (length(dist) == 1) p1   else p1[ind]
      priors$p2[j]    <- if (length(dist) == 1) p2   else p2[ind]
      pr              <- build_prior(dist = priors$dist[j], p1 = priors$p1[j], p2 = priors$p2[j])
      priors$mean[j]  <- pr$mean
      priors$sd[j]    <- pr$sd
      priors$fun[[j]] <- pr$fun
      if (priors$dist[j] == "fixed") ifixed <- c(ifixed, j)
    }
  }


  cache$inotFixed  <- inotFixed  <- if (!is.null(ifixed)) -ifixed else 1:l
  cache$ifixed     <- ifixed
  cache$priors     <- priors
  Phi[1, ]         <- priors$mean
  Phi[, ifixed]    <- matrix(rep(priors$mean[ifixed], nMCMC), nrow = nMCMC, byrow = T)
  cache$sdRates    <- priors$sd


  # data matrices and vectors and scaling
  Xs             <- sim[, 2:(p + q + 1)]
  ys             <- sim[, 1]
  Xf             <- field[, 2:(p + 1), drop = F]
  yf             <- field[, 1]

  meanYs         <- mean(ys)
  sdYs           <- sd(ys)
  expMin         <- apply(Xs[, iexp, drop = F], 2, min)
  expMax         <- apply(Xs[, iexp, drop = F], 2, max)
  calMin         <- apply(Xs[, ical, drop = F], 2, min)
  calMax         <- apply(Xs[, ical, drop = F], 2, max)
  expRange       <- expMax - expMin
  calRange       <- calMax - calMin
  cache$scale    <- list(meanYs = meanYs, sdYs     = sdYs,
                         expMin = expMin, expRange = expRange,
                         calMin = calMin, calRange = calRange)
  Xs[, iexp]     <- matrix(scale(Xs[,iexp], center = expMin, scale = expRange), ncol = length(iexp))
  Xs[, ical]     <- matrix(scale(Xs[,ical], center = calMin, scale = calRange), ncol = length(ical))
  Xf             <- matrix(scale(Xf,        center = expMin, scale = expRange), ncol = length(iexp))
  y              <- as.vector(scale(c(yf, ys), center = meanYs, scale = sdYs))

  TrueKappa      <- matrix(replicate(n, Phi[1, ikappa]), nrow = n, byrow = T)
  Xk             <- cbind(Xf, TrueKappa)

  cache$Xs       <- Xs
  cache$Xk       <- Xk
  cache$Xf       <- Xf
  cache$y        <- y


  # set up first correlation matrices and load them into cache
  # CorKK : (n * n) correlation matrix of augmented Xk's
  # CorKS : (n * m) correlation matrix between Xk's, Xs's
  # CorSK : (m * n) correlation matrix between Xs's, Xk's
  # corSS : (m * m) correlation matrix between Xs's
  # CorFF : (n * n) correlation matrix between Xf's
  cache$CorKK  <- CorKK   <- correlation(Xk, Xk, theta = Phi[1, ithetaS], alpha = Phi[1, ialphaS])
  cache$CorKS  <- CorKS   <- correlation(Xk, Xs, theta = Phi[1, ithetaS], alpha = Phi[1, ialphaS])
  cache$CorSK  <- CorSK   <- t(CorKS)
  cache$CorSS  <- CorSS   <- correlation(Xs, Xs, theta = Phi[1, ithetaS], alpha = Phi[1, ialphaS])
  cache$CorFF  <- CorFF   <- correlation(Xf, Xf, theta = Phi[1, ithetaB], alpha = Phi[1, ialphaB])

  # compute the first log likelihood by forming the augmented covariance matrix
  #   and load both augmented covariance matirx and log likelihood into cache
  cache$Inn <- Inn <- diag(n)
  AugCov    <- rbind(cbind(Phi[1, isigma2S]*CorKK + Phi[1, isigma2B]*CorFF + Phi[1, isigma2E]*Inn,
                           Phi[1, isigma2S]*CorKS),
                     cbind(Phi[1, isigma2S]*CorSK, Phi[1, isigma2S]*CorSS))

  CholCov   <- chol(AugCov)
  InvCov    <- chol2inv(CholCov)
  logDetCov <- 2*sum(log(diag(CholCov)))

  # computes posterior log likelihood of augmented response given the augmented
  #   covariance matrix (its inverse and determinant) and residuals
  logPost        <- double(nMCMC)
  res            <- y - Phi[1, imuB]
  logPrior       <- 0
  for (i in 1:l) {
    fn           <- priors$fun[[i]]
    logPrior     <- logPrior + fn(Phi[1, i])
  }
  logPost[1]     <- logPrior - 0.5*(logDetCov + drop(t(res) %*% InvCov %*% res))
  if (showProgress) {
    report            <- rbind(round(Phi[1, inotFixed], 3), round(cache$priors$sd[inotFixed], 3))
    colnames(report)  <- cache$priors$param[inotFixed]
    row.names(report) <- c("initial value", "initial proposal sd")
    print(report)
  }
  return(list(Phi = Phi, logPost = logPost))
}




