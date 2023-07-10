#' Calibration Object Summary
#'
#' `summary.fbc()` produces a dataframe of summary statistics for all model parameters
#'
#' @param object  `fbc` object
#'
#' @return dataframe where rows represent model parameters and columns represent: mean,
#'         median, mode, lwr50, upr50, lwr80, upr80, sd
#'
#' @example examples/ex_summary.R
#' @export
summary.fbc <- function(object) {
  return(object$estimates)
}
