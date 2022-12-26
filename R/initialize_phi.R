# Initialize the first row of the phi matrix to be used as initialization for
# MCMC algorithm
initialize_phi <- function(sim, k, p, q, theta0, omega0, alpha0, sigma20) {
  phi_init <- double(k)
  phi_init[1 : q] <- theta0
  phi_init[1 + q : q + (p + q)] <- omega0
  phi_init[1 + q + (p + q) : q + (2 * (p + q))] <- alpha0
  phi_init[1 + q + (2 * (p + q)) : q + (2 * (p + q)) + p] <- omega0
  phi_init[1 + q + (2 * (p + q)) + p : q + (2 * (p + q)) + (2 * p)] <- alpha0
  phi_init[1 + q + (2 * (p + q)) + (2 * p) : k - 1] <- sigma20
  R <- correlation(X = sim[, 2: 1 + p + q],
                   scale = phi_init[1 + q : q + (p + q)],
                   smoothness = phi_init[1 + q + (p + q) : q + (2 * (p + q))])
  phi_init[k] <- mu_hat(sim[, 1], R)
  return(phi_init)
}
