sigma2_pr <- function(sigma2, type) {
  if (type=="inverse gamma") {
    return(0.5 * (exp(lambda)/((1+exp(lambda))^2))
           / sqrt(exp(-lambda)/(1+exp(-lambda))))
  }
  return(-1)
}
