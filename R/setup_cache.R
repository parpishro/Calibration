setup_cache <- function(sim, field) {

  # scalers
  m   <- nrow(sim)               # number of simulation runs
  n   <- nrow(field)             # number of field observations
  p   <- ncol(field) - 1         # number of experimental variables
  q   <- ncol(sim) - p - 1       # number of calibration parameters
  d   <- ncol(sim) - 1           # number of all variables for simulation
  k   <- q + (p + q) + (p + q) +  p + p + 1 + 1 + 1 + 1 # total # of parameters

  # index range identifiers
  theta   <- 1:q    # calibration
  omegaS  <- (q+1): (q + (p + q))                       # sim scale
  alphaS  <- (q + (p + q) + 1): (q + (p + q) + (p + q)) # sim smoothness
  omegaB  <- (q + (p + q) + (p + q) + 1): (q + (p + q) + (p + q) + p)#bias scale
  alphaB  <- (q + (p + q) + (p + q) + p + 1): (k - 4)   #  bias smoothness
  sigma2S <- k - 3  # sim variance
  sigma2B <- k - 2  # bias variance
  sigma2E <- k - 1  # error variance
  mu      <- k      # number of total parameters

  # data matrices and vectors
  Xs      <- sim[, 1:d]
  ys      <- sim[, d + 1]
  Xb      <- field[, 1:p, drop = FALSE]
  yf      <- field[, p + 1]
  y       <- (c(ys, yf) - mean(ys)) / sd(ys)

  # parameters (initialize first row of Phi matrix)
  Phi              <- matrix(nrow = Nmcmc, ncol = k)
  Phi[1, theta]    <- apply(Xs[, theta, drop = FALSE], 2, mean)
  Phi[1, omegaS]   <- double(length(omegaS)) + 1
  Phi[1, alphaS]   <- double(length(alphaS)) + 1.8
  Phi[1, omegaB]   <- double(length(omegaB)) + 1
  Phi[1, alphaB]   <- double(length(alphaB)) + 1.8
  Phi[1, sigma2S]  <- 1
  Phi[1, sigma2B]  <- 1
  Phi[1, sigma2E]  <- 1

  CorSS            <- correlation(Xs, Phi[1, omegaS], Phi[1, alphaS])
  Phi[1, mu]    <- mu_hat(CorSS, ys)

  # setting cache environment: TO BE ACCESSED/MODIFIED BY ALL PACKAGE FUNCTIONS
  assign('m', m, envir = cache)    # number of simulation runs
  assign('n', n, envir = cache)    # number of field observations
  assign('p', p, envir = cache)    # number of experimental variables
  assign('q', q, envir = cache)    # number of calibration parameters
  assign('d', d, envir = cache)    # number of all variables for simulation
  assign('k', k, envir = cache)

  assign('theta',   theta ,  envir = cache)     # calibration
  assign('omegaS',  omegaS,  envir = cache)     # sim scale
  assign('alphaS',  alphaS,  envir = cache)     # sim smoothness
  assign('omegaB',  omegaB,  envir = cache)     # bias scale
  assign('alphaB',  alphaB,  envir = cache)     # bias smoothness
  assign('sigma2S', sigma2S, envir = cache)     # sim variance
  assign('sigma2B', sigma2B, envir = cache)     # bias variance
  assign('sigma2E', sigma2E, envir = cache)     # error variance
  assign('mu',      mu,      envir = cache)     # number of total parameters

  assign('Xs', Xs, envir = cache)
  assign('ys', ys, envir = cache)
  assign('Xb', Xb, envir = cache)
  assign('yf', yf, envir = cache)
  assign('y',  y,  envir = cache)

  assign('Phi',   Phi,    envir = cache)
  assign('CorSS', CorSS,  envir = cache)

}
