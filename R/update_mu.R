# compute the mu_hat based on current (last update) correlation function
# parameters

update_mu <- function() {

  f     <- double(length(cache$y)) + 1
  uf    <- backsolve(cache$Chol, f, transpose = TRUE)
  ftrf  <- crossprod(uf)
  uy    <- backsolve(cache$Chol, cache$y, transpose = TRUE)
  ftry  <- crossprod(uf, uy)
  muHat <- solve(ftrf, ftry)
  return(drop(muHat))
}
