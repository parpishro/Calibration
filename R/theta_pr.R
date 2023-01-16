theta_pr <- function(lambda, type = "uniform") {
  if (type=="uniform") {
    return(0.5 * (exp(lambda)/((1+exp(lambda))^2))
           / sqrt(exp(-lambda)/(1+exp(-lambda))))
  }
  return(-1)
}
