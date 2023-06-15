#' Control the optional arguments of `calibrate()`
#'
#' @param ... Any of the optional arguments in `calibrate()` function
#'
#' @export
control <- function(...) {


  args <- list(...)
  for (argName in names(args)) {
    if (argName %in% names(formals(calibrate))) {
      from <- formals(calibrate)[[argName]]
      formals(calibrate)[[argName]] <<- args[[argName]]
      cat("Argument `", argName, "` changed from ", from, " to ", args[[argName]], sep = "")
    }

    else
      warning(paste0("Argument '", argName, "' not found in calibrate() function"))
  }
}
