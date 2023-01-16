#' Conditional (posterior) log likelihood of augmented response vector (z)
#'
#' @param CovD augmented covariance matrix
#'
#' @return posterior log likelihood of augmented response given the augmented
#'            covariance matrix (its inverse and determinant) and residuals
log_lik <- function (chol, res) {
  return (-0.5 * (log(chol$det)) - (res %*% chol$inv %*% res))
}


