#' Calibration Object Print
#'
#' `print.fbc()` prints the prior specification and parameter summary for all model parameters
#'
#' @param x  `fbc` object
#' @param ...      other arguments
#'
#' @example man/examples/ex_print.R
#' @export
print.fbc <- function(x, ...) {
  print(x$estimates)
}
