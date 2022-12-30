proposal <- function(param) {
  if (length(param) == 1) {
    sdSoFar  <- 5
  } else {
    sdSoFar  <- sqrt(var(param))
  }
  last     <- param[length(param)]
  proposed <- 0.95 * rnorm(1, mean = last, sd = 2.38 *sdSoFar)
              + 0.05 * rnorm(1, mean = last, sd = 0.1)
  return(proposed)
}
