# Build a simple calibration model (1 experimental input + 1 calibration input)
cal1    <- calibrate(sim = analytic11S, field = analytic11F,
                     nMCMC = 5, nBurn = 0, thinning = 1,
                     kappaDist = "normal", kappaP1 = 0.5, kappaP2 = 0.25)

# plot default plots the density plot of calibration parameter
plot(cal1)

# plot other model parameters (note that since there are two correlation scale parameters for
# simulator GP, there will be two plots)
par(mfrow = c(1, 2))
plot(cal1, parameter = "thetaS")

# plot the trace of MCMC draws (note that this plot is useful to determine convergence and
# stationarity of the parameter draws. This minimal example dhas only few draws!)
par(mfrow = c(1, 2))
plot(cal1, parameter = "alphaS", type = "trace")


# plot the fitted values of the calibration model. It will plot the fitted values versus all
# experimental inputs, in this example only one!
plot(cal1, type = "fits")

# Build more complex calibration model (2 experimental input + 2 calibration input)
cal2    <- calibrate(sim = Ds2, field = Df2, nMCMC = 20, nBurn = 0, thinning = 1)

# plot default plots the density plot of calibration parameters (two plots for two
# calibration parameters)
plot(cal2)
