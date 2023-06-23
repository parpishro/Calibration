
plot <- function(x) {
  UseMethod("plot")
}


#' Plot Calibration Model
#'
#' Given a parameter, `plot.fbc()` can plot three types of plots:
#'  - "density": plots the prior and posterior density distribution of chosen parameter
#'  - "trace":   plots the trace of the chosen parameter over MCMC runs (afte thinning)
#'  - "fits":    plots the fit of field training data to calibrated model (using "MAP" method)
#'
#' @param object     `fbc` object (output of `calibrate()` function)
#' @param parameter  character string representing the parameter class
#' @param type       type of plotting
#'
#' @export
#'
#' @examples
plot.fbc <- function(object, parameter = "kappa", type = "density") {
  c   <- object$cache
  pr  <- object$priors[parameter][[1]]
  Phi <- object$Phi
  ind <- c$indices[paste0("i", parameter)][[1]]

  if (type == "density") {
    for (i in ind) {
      label <- if (parameter %in% c("sigma2S", "sigma2B", "sigma2E", "muB")) parameter else paste0(parameter, i-ind[1]+1)

      prior_fn <- c$priorFns[[i]]
      prX      <- switch(parameter, kappa = (seq(0.001,1,0.001) * c$scale$calRange[i]) + c$scale$calMin[i], seq(0.002,2,0.002))
      prY      <- switch(parameter, kappa = exp(prior_fn(seq(0.001,1,0.001))), exp(prior_fn(seq(0.002,2,0.002))))
      posX     <- density(Phi[ ,label], cut=0.6)$x
      posY     <- density(Phi[ ,label], cut=0.6)$y
      xlimit   <- c(min(c(prX, posX)), max(c(prX, posX))) * c(0.95, 1.05)
      ylimit   <- c(0, max(c(prY, posY))) * 1.2
      pMode    <- object$estimates[label, 'mode']

      plot(posX, posY, main="Density Plot", xlab=label, ylab= "density", type="l", lwd=2, col="blue", xlim=xlimit, ylim=ylimit)
      lines(prX, prY,  type="l")
      abline(v=pMode, col="lightblue", lty=2)
      legend("top", c("Posterior", "Prior", "Mode"), lty=c(1,1,2), text.font = 4, col=c("blue", "black", "lightblue"),  cex = 0.4, horiz=T)
    }
  } else if (type == "trace") {
    for (i in ind) {
      label <- if (parameter %in% c("sigma2S", "sigma2B", "sigma2E", "muB")) parameter else paste0(parameter, i-ind[1]+1)
      plot(Phi[, i], type="l", xlab="index", ylab=label, main="Trace Plot")
    }
  } else if (type == "fits") {
    X       <- matrix(scale(unique(c$Xf), center=-(c$scale$expMin/c$scale$expRange), scale=1/c$scale$expRange), ncol = length(ind), byrow = T)
    preds   <- predict.fbc(object = object, newdata = X, type = "MAP")
    actuals <- (c$y[1:c$n]*c$scale$sdYs)+c$scale$meanYs
    fits    <- preds$pred
    lower95 <- fits-(2*preds$se)
    upper95 <- fits+(2*preds$se)

    for (i in 1:ncol(c$Xf)) {

      ylimit <- c(min(c(lower95, actuals)), max(c(upper95, actuals))) * c(0.8, 1.2)
      xf <- (c$Xf[, i] * c$scale$expRange[i]) + c$scale$expMin[i]
      plot(xf,      actuals, cex = 0.65, pch = 1,  col = "black", main="Fits Plot", xlab=paste0("x", i), ylab="field response", ylim=ylimit)
      lines(X[, i], fits,    cex = 0.65, pch = 15, col = "blue",      type="o")
      lines(X[, i], upper95, cex = 0.65, pch = 0,  col = "lightblue", type="o", lty=2)
      lines(X[, i], lower95, cex = 0.65, pch = 0,  col = "lightblue", type="o", lty=2)
      legend("top", legend = c("Field Response", "Prediction", "95% P.I."), lty=c(0,1,2), text.font = 4,
             col = c("black", "blue", "lightblue"), pch = c(1, 15, 0), cex = 0.5, horiz=T)
    }
  } else {
    stop("Invalid plot type!")
  }
}
