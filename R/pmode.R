#' Compute Mode of a Vector
#'
#' @param x       numeric vector
#' @param breaks  integer representing number of bins to compute mode
#'
#' @return double represent estimated mode of
#' @export
#'
#' @examples
pmode <- function(x, breaks = 20) {
  bins   <- seq(min(x), max(x), length.out = breaks+1)
  counts <- double(breaks)
  for (i in 2:(beaks+1))
    counts[i-1] <- sum(x >= bins[i-1] & x < bins[i])
  return(mean(bins[which.max(counts):(which.max(counts)+1)]))
}
