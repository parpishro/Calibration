#' Estimate Mode from a Sample
#'
#' `pmode` bins sample observations into equally spaced intervals and counts number of
#' observations in each bin. The bin with highest number of observations represents the
#' mode and mean of observations in that bin is returned as representative mode of the
#' vector.
#'
#' @param x       numeric vector
#' @param breaks  integer representing number of bins to compute mode
#'
#' @return double represent estimated mode of
#' @export
#'
#' @examples
pmode <- function(x, breaks = NA) {
  if (length(x) < 6)
    return(mean(x))
  if (is.na(breaks))
    breaks <- floor(length(x)/5)
  bounds <- seq(min(x), max(x), length.out = breaks+1)
  counts <- double(breaks)
  bins   <- list()
  for (i in 1:breaks) {
    bins[[i]] <- x[x >= bounds[i] & x < bounds[i + 1]]
    counts[i] <- length(bins[[i]])
  }
  counts[breaks] <- counts[breaks] + 1
  return(mean(bins[[which.max(counts)]]))
}
