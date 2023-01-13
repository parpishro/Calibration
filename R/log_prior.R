log_prior <- function(theta, lambda, gamma, sigma2) {
  return(sum(log(theta_pr(theta))) +
         sum(log(lambda_pr(lambda))) +
         sum(log(gamma_pr(gamma))) +
         sum(log(sigma2_pr(sigma2))))
}
