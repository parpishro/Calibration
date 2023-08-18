#' `fbc` Class Constructor
#'
#' Creates an object of class `fbc` from a given list with correct elements
#'
#' @param x     named list of length 7 containing "Phi", "estimates", "logPost", "priors",
#'              "acceptance", "vars", "cache" as fields
#'
#' @return object of class `fbc` object with mentioned elements`
#' @export
fbc <- function(x) {
  stopifnot(identical(names(x),
                      c("Phi", "estimates", "logPost", "priors", "acceptance", "vars", "data",
                        "scale", "indices")))
  return(structure(x, class = "fbc"))
}
