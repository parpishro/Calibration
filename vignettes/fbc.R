## ---- include = FALSE-----------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ">"
)
options(width = 100)

## ----install FBC, eval=FALSE----------------------------------------------------------------------
#  install.packages(FBC)

## ----install development version of FBC, eval=FALSE-----------------------------------------------
#  devtools::install_github("parpishro/FBC")

## ----load FBC-------------------------------------------------------------------------------------
library(FBC)

## ----ball data matrices, echo=TRUE----------------------------------------------------------------
head(ballField, 3)
dim(ballField)
head(ballSim, 3)
dim(ballSim)

## ----ball plot html, echo=FALSE, fig.height=5, fig.width=7----------------------------------------
plot(ballField[, 2], ballField[, 1], cex = 0.65, pch = 19, col = "red", ylim = c(0, 1.5), 
     xlab = "Height (m)", ylab = "Time (s)")
points(ballSim[, 2], ballSim[, 1], cex = 0.65, pch = 19, col = "blue")
legend("topleft", legend = c("Simulation", "Field"), 
       col = c("blue", "red"), pch = 16, cex = 0.8)

## ----calibrate, eval=FALSE------------------------------------------------------------------------
#  calMod <- calibrate(sim = ballSim, field = ballField,                 # Data
#                      nMCMC = 11000, nBurn = 1000, thinning = 50,       # MCMC
#                      kappaDist = "beta", kappaP1 = 1.1, kappaP2 = 1.1, # Priors
#                      hypers = set_hyperPriors(),
#                      showProgress = FALSE)

## ----load calMod, include=FALSE-------------------------------------------------------------------
calMod_path <- system.file("calMod.rda", package = "FBC", mustWork = TRUE)
load(calMod_path)

## ----calMod---------------------------------------------------------------------------------------
names(calMod)

## ----Phi columns----------------------------------------------------------------------------------
head(calMod$Phi, 3)

## ----predict--------------------------------------------------------------------------------------
predsMAP   <- predict(object = calMod, newdata = matrix(c(2.2, 2.4), ncol = 1), method = "MAP")
predsBayes <- predict(object = calMod, newdata = matrix(c(2.2, 2.4), ncol = 1), method = "Bayesian")

## -------------------------------------------------------------------------------------------------
predsMAP  
predsBayes

## ----set_hyperPriors------------------------------------------------------------------------------
priors <- set_hyperPriors(thetaSDist = "beta", thetaBP2 = 6)

## ----plot density---------------------------------------------------------------------------------
# Note that there are two correlation scale parameters in simulator GP and there will be two plots
plot(calMod, parameter = "thetaS")

## ----plot trace-----------------------------------------------------------------------------------
# Note that there are two correlation smoothness parameters in simulator GP and there will be 
# two plots
plot(calMod, parameter = "alphaS", type = "trace")

## ----plot fits------------------------------------------------------------------------------------
# Plots the fitted values versus all experimental inputs along with actual values in separate plots
plot(calMod, type = "fits", xlab = "height")

## ----summary and print----------------------------------------------------------------------------
calModSum <- summary(calMod) 
print(calMod)

## ----analytic11-----------------------------------------------------------------------------------
# priors    <- set_hyperPriors(thetaSDist  = "logbeta",      thetaSP1  = 5,  thetaSP2 = 5,
#                              alphaSDist  = "fixed",        alphaSP1  = 2,
#                              thetaBDist  = "logbeta",      thetaBP1  = 5,  thetaBP2 = 5,
#                              alphaBDist  = "fixed",        alphaBP1  = 2,
#                              sigma2SDist = "inversegamma", sigma2SP1 = 5,  sigma2SP2 = 0.2,
#                              sigma2BDist = "inversegamma", sigma2BP1 = 12, sigma2BP2 = 0.5,
#                              sigma2EDist = "inversegamma", sigma2EP1 = 5,  sigma2EP2 = 100,
#                              muBDist     = "fixed",        muBP1     = 0)
# analytic1 <- calibrate(analytic11S, analytic11F, nMCMC = 15000, nBurn = 5000, 
#                        kappaDist = "normal", kappaP1 = 0.5, kappaP2 = 0.25, hypers = priors)

