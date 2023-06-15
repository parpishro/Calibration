#' predict
#'
#' @param object     class fbc
#' @param newdata a vector of new data
#'
#' @return predictions
#'
#' @export
predict.fbc <- function(object, newdata, type) {
  cache     <- object@cache
  Phi       <- object@Phi
  InvCovRes <- cache$InvCovRes
  expMin    <- cache$scale$expMin
  expRange  <- cache$scale$expRange
  calMin    <- cache$scale$calMin
  calRange  <- cache$scale$calRange
  yMean     <- cache$scale$meanYs
  ySd       <- cache$scale$sdYs
  Ystar     <- matrix(0, nrow = nrow(Phi), ncol = nrow(newdata))
  Xstar     <<- matrix(scale(newdata,
                               center=expMin,
                               scale=expRange), ncol=length(expMin))
  Phi[,cache$ikappa] <- matrix(scale(Phi[,cache$ikappa],
                                     center=calMin,
                                     scale=calRange), ncol=length(cache$ikappa))

  if (type == "MAP") {
    iMAP        <- which.max(object@logPost)
    phiMAP      <- as.numeric(Phi[iMAP, ,drop=T])
    Xk          <- cbind(cache$Xf, replicate(cache$n, phiMAP[cache$ikappa]))
    Xkstar      <- cbind(Xstar,    replicate(cache$n, phiMAP[cache$ikappa]))
    CorSN       <- correlation(cache$Xs, Xkstar, theta=phiMAP[cache$ithetaS], alpha=phiMAP[cache$ialphaS])
    CorKN       <- correlation(Xk,        Xkstar, theta=phiMAP[cache$ithetaS], alpha=phiMAP[cache$ialphaS])
    CorFN       <- correlation(cache$Xf, Xstar,  theta=phiMAP[cache$ithetaB], alpha=phiMAP[cache$ialphaB])
    CovN        <- rbind(phiMAP[cache$isigma2S]*CorKN + phiMAP[cache$isigma2B]*CorFN, phiMAP[cache$isigma2S]*CorSN)
    Ystar       <- phiMAP[cache$imu] + (t(CovN) %*% InvCovRes[iMAP, ])
  } else if (type == "Bayesian") {
    Inn       <- diag(cache$n)
    cat("Statrting prediction draws ... \n")
    Phi[,cache$ikappa] <- matrix(scale(Phi[,cache$ikappa],
                                       center=calMin,
                                       scale=calRange), ncol=length(cache$ikappa))
    for (i in 1:nrow(Phi)) {
      cache$Xk    <- cbind(cache$Xf, replicate(cache$n, Phi[i, cache$ikappa]))
      Xkstar      <- cbind(Xstar,    replicate(cache$n, Phi[i, cache$ikappa]))
      CorSN       <- correlation(cache$Xs, Xkstar, theta=Phi[i,cache$ithetaS], alpha=Phi[i,cache$ialphaS])
      CorKN       <- correlation(cache$Xk, Xkstar, theta=Phi[i,cache$ithetaS], alpha=Phi[i,cache$ialphaS])
      CorFN       <- correlation(cache$Xf, Xstar,  theta=Phi[i,cache$ithetaB], alpha=Phi[i,cache$ialphaB])
      CovN        <- rbind(Phi[i,cache$isigma2S]*CorKN + Phi[i,cache$isigma2B]*CorFN, Phi[i, cache$isigma2S]*CorSN)
      Ystar[i,]   <- Phi[i,cache$imu] + (t(CovN) %*% InvCovRes[i, ])
      if ((i %% (nrow(Phi)/10)) == 0)
        cat("Completed ", 100*(i/nrow(Phi)), " % of prediction draws ... \n")
      }
  } else
    stop("Unknown prediction type!")
  return((Ystar * ySd) + yMean)
}


