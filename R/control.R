#' Control the optional arguments of `calibrate()`
#'
#' @param ... Any of the optional arguments in `calibrate()` function
#'
#' @export
control <- function(...) {
  args <- list(...)
  for (argName in names(args)) {
    if (argName %in% names(formals(calibrate))) {
      formals(calibrate)[[argName]] <- args[[argName]]
      print(c(formals(calibrate)[[argName]], args[[argName]]))
    }

    else
      warning(paste0("Argument '", argName, "' not found in calibrate() function"))
  }
}
