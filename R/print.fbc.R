#' Calibration Object Print
#'
#' `print.fbc()` prints the prior specification and parameter summary for all model parameters
#'
#' @param x  `fbc` object
#'
#' @example examples/ex_print.R
#' @export
print.fbc <- function(x) {
  print(x$estimates)
}
