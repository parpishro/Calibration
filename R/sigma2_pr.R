sigma2_pr <- function(sigma2, type = "inverse gamma") {
  if (type=="inverse gamma") {
    return(0.5 * (exp(sigma2)/((1+exp(sigma2))^2))
           / sqrt(exp(-sigma2)/(1+exp(-sigma2))))
  }
  return(-1)
}
