# compute the mu_hat based on current (last update) correlation function
# parameters

update_mu <- function() {

  if (!is.null(cache$InvCov)) {
    ones  <- matrix(1, nrow = cache$n + cache$m)
    u     <- crossprod(ones, cache$InvCov)
    term1 <- drop(tcrossprod(u, t(ones)))
    term2 <- drop(tcrossprod(u, t(cache$y)))
    muHat <- term2 / term1
  } else if (exists('muHat', envir = cache))
    muHat <- cache$muHat
  return(muHat)
}

