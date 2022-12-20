#' Conditional log likelihood of augmented response vector (d)
#'
#' @param CovD_inv inverse of augmented covariance matrix
#' @param CovD_det determinant of augmented covariance matrix
#' @param d augmented response vector which is (yT, zT)T
#' @param mu_y mean of simulator response vector (y)
#'
#' @return log likelihood of augmented response given the parameters phi
#' @export
#'
#' @examples
ll_d <- function (CovD_inv, CovD_det, d, mu_y) {
  r <- d - mu_y
  return (-0.5 * (log(CovD_det) - (r %*% CovD_inv %*% r)))
}
