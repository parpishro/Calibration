#' Conditional log likelihood of augmented response vector (d)
#'
#' @param CovD_inv inverse of augmented covariance matrix
#' @param CovD_det determinant of augmented covariance matrix
#' @param d augmented response vector which is (yT, zT)T
#' @param mu_y mean of simulator response vector (y)
#'
#' @return log likelihood of augmented response given the parameters phi
log_lik <- function (logLik, phi, chnaged, env) {

  env0             <- environment()
  parent.env(env0) <- env

  covD             <- update_cov(phi, changed, env)



  r <- d - mu_hat
  return (-0.5 * (log(CovD_det) - (r %*% CovD_inv %*% r)))
}
