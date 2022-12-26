# Calibrates the simulator using both simulator and field data and returns the
# distribution of calibration parameters and hyperparameters of the model


# REQUIRE: simulation data (sim) must be a matrix with rows representing
#             simulator observations and first column being the univariate
#             response variable, followed by experimental input columns,
#             followed by calibration input columns.
#          field data (field) must be a matrix with rows representing
#             field observations and first column being the univariate response
#             variable, followed by experimental input columns.
#          prior type for calibration parameters (theta_pr) must be either
#             "weak" (default) or "strong".
#          prior for scale of correlation functions (omega_pr) must be "logbeta"
#          prior for smoothness correlation function (alpha_pr) must be
#             "logistic"
#          prior for uncertainty (sigma2_pr) must be "inverse gamma"
# EFFECT: given simulation data, field data, priors for calibration parameters,
#             correlation functions hyperparameters, uncertainty hyperparameters
#             , and MCMC specifications, build a KOH model and run MH MCMC to
#             find the posterior distribution of parameters, hyperparameters,
#             and calibrated response.

calibrate <- function(sim, field,
                      Nmcmc = 10000, nBurn = 500, thining = 100,
                      theta_pr = "weak", omega_pr = "logbeta",
                      alpha_pr = "logistic", sigma2_pr = "inverse gamma",
                      theta0, omega0, alpha0, sigma20) {
  clb <- mcmc(sim, field,
              Nmcmc, nBurn, thining,
              theta_pr, omega_pr, alpha_pr, sigma2_pr,
              theta0, omega0, alpha0, sigma20)
  paramVar <- apply(clb$mean,2,var) + apply(clb$var,2,mean) # law of total variance
  paramMean <- apply(clb$mean,2,mean)
  results <- list(mean = paramMean, var = paramVar)
  return(results)
}


runBMCMC<-function(nmcmc,burn,thin,x,y,xtest1, lambda.ini, lambda.w.ini, gamma.ini, gamma.w.ini){

  m<-list(pred.y=pred.y, pred.var=pred.var, reasonable.lambda=xxy$reasonable.lambda, accept.lambda=xxy$accept.lambda,
          mcmc.ma.lambda=xxy$mcmc.ma.lambda,
          accept.p=xxy$accept.gamma,
          mcmc.ma.gamma=xxy$mcmc.ma.gamma,
          mcmc.ma.p=1+1/(1+exp(-(xxy$mcmc.ma.gamma))) )
  return(m)
}
