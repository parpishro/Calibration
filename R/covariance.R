#' Compute the augmented covariance matrix
#'
#' @param sim simulation data
#' @param field field data
#' @param phi model parameters/hyperprameters
#'
#' @return covariance matrix of augmented data matrix
#' @export
#'
#' @examples
covariance <- function(sim, field, phi) {

  m        <- nrow(sim)
  pq        <- ncol(sim)
  n       <- nrow(field)
  p        <- ncol(field)

  Xs     <- sim[, 2:ks]
  Ys      <- sim[, 1]
  Xb       <- field[, 2:p]
  Xf    <- cbind(Xb, matrix(replicate(phi[1:(pq - p)]), n))
  Yf    <- field[, 1]


  s2s     <- phi[len(phi) - 3] # marginal variance of simulator
  s2b     <- phi[len(phi) - 2] # marginal variance of bias correction
  s2e     <- phi[len(phi) - 1] # marginal variance of noise (measurement error)



  # (n * n) correlation matrix between augmented Xf's
  CorFF     <- correlation(Xfield)

  CorFS     <- correlation(Xfield, Xsim) # (n * m) correlation matrix between Xf's, Xs's
  CorSF     <- t(CorFS)          # (m * n) correlation matrix between Xs's, Xf's
  CorSS     <- correlation(Xsim)         # (m * m) correlation matrix between Xs's
  CorB   <- correlation()                # (n * n) correlation matrix between Xb's



  Inn       <- diag(n)

  CovD <- cbind(rbind((s2s * CorFF) + (s2b * XX) + (s2e * Inn), CorFS),
                rbind((s2s * CorSF), (s2s * CorSS)))

  return(CovD)
}


