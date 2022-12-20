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
func <- function(RS_XsT, RS_XsT_XfTheta, RS_XfTheta, RB_Xf, vs, vb, ve) {
  n <- nrow(RS_XsT)
  m <- nrow(RS_XfTheta)
  I_nn  <- diag(n)
  zero_mm <- matrix(0, m, m)
  zero_mn <- matrix(0, m, n)
  zero_nm <- matrix(0, n, m)
  CS <- vs * rbind(cbind(RS_XsT, RS_XsT_XfTheta),
                   cbind(t(RS_XsT_XfTheta), RS_XfTheta))
  CB <- vb * rbind(cbind(zero_mm, zero_mn), cbind(zero_nm, RB_Xf))
  CE <- ve * rbind(cbind(one_nn, zero_nm), cbind(zero_mn, zero_mm))
  return(CovSim + CovDelta + CovEps)
}
