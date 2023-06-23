#' fbc Class
#'
#' @slot estimates data.frame.
#' @slot Phi data.frame.
#' @slot logPost numeric.
#' @slot acceptance numeric.
#' @slot vars character.
#' @slot cache environment.
#'
#' @export
setClass("fbc",
         slots = list(Phi           = "data.frame",
                      estimates     = "data.frame",
                      logPost       = "numeric",
                      priors        = "list",
                      acceptance    = "numeric",
                      vars          = "character",
                      cache         = "environment"))

setMethod("predict",
          signature(object = "fbc"),
          function(object, ...) {
            args <- list(...)
            if (length(args) == 1)
              return(predict.fbc(object, args[[1]]))
            else
              return(predict.fbc(object, args[[1]], args[[2]]))
          }
)

setMethod("plot",
          signature(x = "fbc"),
          function(x, y, ...) {
            args <- list(...)
            plot.fbc(x, y, args[[1]])
          }
)


setMethod("plot",
          signature(x = "fbc"),
          function(x, y, ...) {
            args <- list(...)
            plot.fbc(x, y, args[[1]])
          }
)

setMethod("summary",
          signature(object = "fbc"),
          function(object, ...) {
            summary.fbc(object)
          }
)

setMethod("print",
          signature(x = "fbc"),
          function(x, ...) {
            print.fbc(x)
          }
)


