# build a simple calibration model
cal1    <- calibrate(sim = Ds1, field = Df1,
                     Nmcmc = 5, nBurn = 0, thinning = 1,
                     kappa = list(dist = "gaussian", init = 0.5, p1 = 0.5, p2 = 0.25))
# predict the fitted values (Bayesian method)
preds <- predict(cal1, newdata = matrix(Df1[,2], ncol = 1))

# view the fitted value and their MCMC-based standard error
preds$pred
preds$se

# build a more complex calibration model and predict the fitted values
cal2   <- calibrate(sim = Ds2, field = Df2, Nmcmc = 10, nBurn = 0, thinning = 1)
predsB <- predict(cal2, newdata = Df2[, 2:3], method = "Bayesian")
predsB$pred
predsB$se

# and predict using MAP method (much faster)
predsM <- predict(cal2, newdata = Df2[, 2:3], method = "MAP")
predsM$pred
predsM$se
