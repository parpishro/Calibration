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
         slots = list(estimates     = "data.frame",
                      Phi           = "data.frame",
                      logPost       = "numeric",
                      acceptance    = "numeric",
                      vars          = "character",
                      cache         = "environment"))

setMethod("predict",
          signature(object = "fbc"),
          function(object, ...) {
            args <- list(...)
            if(length(args) == 1)
              args <- c(args, "MAP")
            return(predict.fbc(object, args[[1]], args[[2]]))
          }
)


