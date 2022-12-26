# compute the mu_hat based on current (last update) correlation function
# parameters

mu_hat <- function(y, R) {
  F <- double(length(y)) + 1
  invR <- chol2inv(chol(R))
  return((crossprod(F, invR) %*% y) / chol2inv(chol(crossprod(F, invR) %*% F)))
}
