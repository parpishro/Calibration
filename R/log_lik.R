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
log_lik <- function (CovD_inv, CovD_det, d, mu_hat, chnaged = NULL) {

  if (changed == "calibration") {

  } else if (changed == "omega sim") {

  } else if (changed == "alpha sim") {

  } else if (changed == "omega bias") {

  } else if (changed == "alpha bias") {

  } else if (changed == "variance sim") {

  } else if (changed == "variance bias") {

  } else if (changed == "variance sim") {

  } else if (is.null(changed)) {

  } else stop("invalid parameter type!")


  r <- d - mu_hat
  return (-0.5 * (log(CovD_det) - (r %*% CovD_inv %*% r)))
}
