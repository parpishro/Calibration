# Initialize the first row of the phi matrix to be used as initialization for
# MCMC algorithm
initialize_phi <- function(env) {

  env0             <- environment()
  parent.env(env0) <- env

  # TODO: generalize prior mean computation
  phiInit          <- double(k)
  phiInit[calib]   <- apply(Xs[, calib, drop = FALSE], 2, mean)
  phiInit[scaleS]  <- double(length(scaleS)) + 1
  phiInit[smoothS] <- double(length(smoothS)) + 1.8
  phiInit[scaleB]  <- double(length(scaleB)) + 1
  phiInit[smoothB] <- double(length(smoothB)) + 1.8
  phiInit[sig2S]   <- 1
  phiInit[sig2B]   <- 1
  phiInit[sig2E]   <- 1

  CorSim           <- corelation(Xs, Xs, phiInit[scaleS], phiInit[smoothS])
  phiInit[muHat]  <- mu_hat(CorSim, ys)

  return(phiInit)
}
