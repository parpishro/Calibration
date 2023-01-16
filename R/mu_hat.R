# compute the mu_hat based on current (last update) correlation function
# parameters

mu_hat <- function(R, y) {

  f    <- double(length(y)) + 1
  invR <- chol2inv(chol(R, pivot = TRUE))
  return((crossprod(f, invR) %*% y) / chol2inv(chol(crossprod(f, invR) %*% f, pivot = TRUE)))
}
