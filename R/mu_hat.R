# compute the mu_hat based on current (last update) correlation function
# parameters

mu_hat <- function(R, y) {

  f <- double(m) + 1
  invR <- chol2inv(chol(R))
  return((crossprod(f, invR) %*% y) / chol2inv(chol(crossprod(f, invR) %*% F)))
}
