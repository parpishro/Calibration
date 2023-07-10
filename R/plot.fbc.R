#' Plot Calibration Model
#'
#' Given a parameter, `plot.fbc()` can plot three types of plots:
#'  - "density": plots the prior and posterior density distribution of chosen parameter
#'  - "trace":   plots the trace of the chosen parameter over MCMC runs (afte thinning)
#'  - "fits":    plots the fit of field training data to calibrated model (using "MAP" method)
#'
#' @param x         `fbc` object (output of `calibrate()` function)
#' @param parameter  character string representing the parameter class
#' @param type       type of plotting
#' @param xlab       character string (or vector) representing the x axis labels for "fits" plot
#' @param ...        other arguments
#'
#' @export
#' @import graphics
#' @importFrom stats density
#'
#' @example man/examples/ex_plot.R
#' @export
plot.fbc <- function(x, parameter = "kappa", type = "density", xlab = NULL, ...) {

  stopifnot(parameter %in% c("kappa", "thetaS", "alphaS", "thetaB", "alphaB", "sigma2S",
                             "sigma2B", "sigma2E"))
  stopifnot(type %in% c("density", "trace", "fits"))

  obj       <- x
  Phi       <- obj$Phi
  estimates <- obj$estimates
  priorFns  <- obj$priorFns
  indices   <- obj$indices
  scale     <- obj$scale
  ind       <- indices[paste0("i", parameter)][[1]]

  if (type == "density") {
    for (i in ind) {
      if (parameter %in% c("sigma2S", "sigma2B", "sigma2E", "muB"))
        label    <- parameter
      else
        label    <- paste0(parameter, i - ind[1] + 1)

      prior_fn <- priorFns[[i]]
      prX      <- switch(parameter,
                         kappa  = (seq(0,1,0.001) * scale$calRange[i]) + scale$calMin[i],
                         thetaS = seq(0,5,0.001),
                         thetaB = seq(0,5,0.001),
                         alphaS = seq(1,2,0.001),
                         alphaB = seq(1,2,0.001),
                         muB    = seq(-1,1,0.001),
                         seq(0,5,0.001))
      prY      <- switch(parameter,
                         kappa  = exp(prior_fn(seq(0,1,0.001)))/scale$calRange[i],
                         thetaS = exp(prior_fn(seq(0,5,0.001))),
                         thetaB = exp(prior_fn(seq(0,5,0.001))),
                         alphaS = exp(prior_fn(seq(1,2,0.001))),
                         alphaB = exp(prior_fn(seq(1,2,0.001))),
                         muB    = exp(prior_fn(seq(-1,1,0.001))),
                         exp(prior_fn(seq(0,5,0.001))))
      prY[1]   <- 0
      posX     <- density(Phi[ ,label], from = min(prX), to = max(prX))$x
      posY     <- density(Phi[ ,label], from = min(prX), to = max(prX))$y
      xlimit   <- c(min(c(prX, posX)), max(c(prX, posX))) * c(0.95, 1.05)
      ylimit   <- c(0, min(max(c(prY, posY)), 100)) * 1.2
      pMode    <- estimates[label, 'mode']
      plot(posX, posY,
           main = "Density Plot", xlab = label, ylab = "density", xlim = xlimit, ylim = ylimit,
           type = "l", lwd = 2, col = "blue")
      lines(prX, prY,  type = "l")
      abline(v = pMode, col = "lightblue", lty = 2)
      legend("top", c("Posterior", "Prior", "Mode"),
             lty = c(1,1,2), col = c("blue", "black", "lightblue"), cex = 0.4, horiz = T)
    }
  } else if (type == "trace") {
    for (i in ind) {
      if (parameter %in% c("sigma2S", "sigma2B", "sigma2E", "muB"))
        label <- parameter
      else
        label <- paste0(parameter, i - ind[1] + 1)
      plot(Phi[, i], type = "l", xlab = "index", ylab = label, main = "Trace Plot")
    }
  } else if (type == "fits") {
    X       <- matrix(scale(unique(obj$data$Xf),
                            center = -(scale$expMin/scale$expRange),
                            scale = 1/scale$expRange), ncol = ncol(obj$data$Xf),
                      byrow = T)
    preds   <- predict.fbc(object = obj, newdata = X, method = "MAP")
    actuals <- (obj$data$y[1:indices$n] * scale$sdYs) + scale$meanYs
    fits    <- preds$pred
    lower95 <- fits - (2*preds$se)
    upper95 <- fits + (2*preds$se)

    xlab    <- if (is.null(xlab)) paste0("x", 1:ncol(obj$data$Xf))
    stopifnot(length(xlab) == ncol(obj$data$Xf))

    for (i in 1:ncol(obj$data$Xf)) {

      ylimit <- c(min(c(lower95, actuals)), max(c(upper95, actuals))) * c(0.8, 1.2)
      xf     <- (obj$data$Xf[, i] * scale$expRange[i]) + scale$expMin[i]

      plot(xf, actuals,
           main = "Fits Plot",
           xlab = xlab[i], ylab = "field response", ylim = ylimit,
           cex = 0.65, pch = 1,  col = "black")
      lines(X[, i], fits,    cex = 0.65, pch = 15, col = "blue",      type = "o")
      lines(X[, i], upper95, cex = 0.65, pch = 0,  col = "lightblue", type = "o", lty = 2)
      lines(X[, i], lower95, cex = 0.65, pch = 0,  col = "lightblue", type = "o", lty = 2)
      legend("top", legend = c("Field Response", "Prediction", "95% P.I."),
             col = c("black", "blue", "lightblue"),
             lty = c(0,1,2), pch = c(1, 15, 0), cex = 0.4, horiz = T)
    }
  } else {
    stop("Invalid plot type!")
  }
}
