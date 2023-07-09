#' Compute Inverse of Augmented Covariance Matrix
#'
#' @param phi  numeric vector of all model parameters
#' @param Xf   numeric matrix of training field input matrix, where each row is a field
#'             observation and each column an experimental inputs
#' @param Xs   numeric matrix of training simulation input matrix, where each row is a
#'             simulator run and each column either an experimental or calibration inputs
#' @param Xk   numeric matrix of augmented training field input matrix, where each row is
#'             a field observation and each column either an experimental input or true
#'             calibration input. Note that true calibration inputs are unknown and will
#'             be replaced by MCMC draws.
#' @param ind  useful indices to access elements of `phi`
#'
#' @return numeric matrix representing the inverse of augmented covariance matrix
#' @noRd
compute_covs <- function(phi, Xf, Xs, Xk, ind) {
  CorKK  <- correlation(Xk,     theta=phi[ind$ithetaS], alpha=phi[ind$ialphaS])
  CorKS  <- correlation(Xk, Xs, theta=phi[ind$ithetaS], alpha=phi[ind$ialphaS])
  CorSK  <- t(CorKS)
  CorSS  <- correlation(Xs,     theta=phi[ind$ithetaS], alpha=phi[ind$ialphaS])
  CorFF  <- correlation(Xf,     theta=phi[ind$ithetaB], alpha=phi[ind$ialphaB])
  AugCov <- rbind(cbind(phi[ind$isigma2S]*CorKK + phi[ind$isigma2B]*CorFF + phi[ind$isigma2E]*diag(ind$n),
                        phi[ind$isigma2S]*CorKS),
                  cbind(phi[ind$isigma2S]*CorSK,
                        phi[ind$isigma2S]*CorSS + 0.000001*diag(ind$m)))

  InvCov <- chol2inv(chol(AugCov))
  return(InvCov)
}
