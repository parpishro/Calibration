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


  m      <- nrow(sim)
  pq     <- ncol(sim)
  n      <- nrow(field)
  p      <- ncol(field)

  Xs     <- sim[, 2:pq]
  Ys     <- sim[, 1]
  Xb     <- field[, 2:p]
  Xf     <- cbind(Xb, matrix(replicate(phi[1:(pq - p)], n), nrow = n))
  Yf     <- field[, 1]


  s2s    <- phi[length(phi) - 3] # marginal variance of simulator
  s2b    <- phi[length(phi) - 2] # marginal variance of bias correction
  s2e    <- phi[length(phi) - 1] # marginal variance of noise (measurement error)



  # (n * n) correlation matrix between augmented Xf's
  CorFF  <- correlation(Xf)

  # (n * m) correlation matrix between Xf's, Xs's
  CorFS  <- correlation(Xf, Xs)

  # (m * n) correlation matrix between Xs's, Xf's
  CorSF  <- t(CorFS)

  # (m * m) correlation matrix between Xs's
  CorSS  <- correlation(Xsim)

  # (n * n) correlation matrix between Xb's
  CorB   <- correlation()

  Inn    <- diag(n)

  CovD   <- cbind(rbind((s2s * CorFF) + (s2b * XX) + (s2e * Inn), CorFS),
                rbind((s2s * CorSF), (s2s * CorSS)))

  return(CovD)
}


