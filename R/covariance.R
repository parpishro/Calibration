#' Compute the augmented covariance matrix
#'
#' @param RS_XsT (m * m) correlation matrix between (xs, t)s
#' @param RS_XsT_XfTheta (m * n) correlation matrix between (xs, t) & (xp, theta)
#' @param RS_XfTheta (n * n) correlation matrix between (xp, theta)s
#' @param RB_Xf (m * m) correlation matrix between xp terms
#' @param vs marginal variance of simulator
#' @param vb marginal variance of bias correction
#' @param ve marginal variance of noise
#'
#' @return covariance matrix of augmented data matrix
#' @export
#'
#' @examples
covariance <- function(sim, field, phi) {

  ns        <- nrow(sim)
  ks        <- ncol(sim)
  nf        <- nrow(field)
  kf        <- ncol(field)

  Xsim      <- sim[, 2:ks]
  Ysim      <- sim[, 1]
  XX        <- field[, 2:kf]
  Xfield    <- cbind(XX, matrix(replicate(phi[1:(ks - kf)]), nf))
  Yfield    <- field[, 1]


  CorFF     <- cor(Xfield)
  CorFS     <- cor(Xfield, Xsim)
  CorSF     <- t(CorFS)
  CorSS     <- cor(Xsim)
  CorBias   <- cor()

  sig2S     <- phi[len(phi) - 3]
  sig2B     <- phi[len(phi) - 2]
  sig2E     <- phi[len(phi) - 1]

  Iff       <- diag(nf)

  CovD <- cbind(rbind((sig2S * CorFF) + (sig2B * XX) + (sig2E * Iff), CorFS),
                rbind((sig2S * CorSF), (sig2S * CorSS)))

  return(CovD)
}
