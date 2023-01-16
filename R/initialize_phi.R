# Initialize the first row of the phi matrix to be used as initialization for
# MCMC algorithm
initialize_phi <- function(env) {

  # TODO: generalize prior mean computation
  phiInit          <- double(k)
  phiInit[calib]   <- apply(Xs[, calib], 2, mean)
  phiInit[scaleS]  <- double(scaleS) + 1
  phiInit[smoothS] <- double(smoothS) + 1.8
  phiInit[scaleB]  <- double(scaleB) + 1
  phiInit[smoothB] <- double(smoothB) + 1.8
  phiInit[sig2S]   <- 1
  phiInit[sig2B]   <- 1
  phiInit[sig2E]   <- 1

  CorSim           <- correlation(X = Xs, phiInit[scaleS], phiInit[smoothS])
  phiInit[muHat]  <- mu_hat(ys, CorSim)

  return(phiInit)
}
