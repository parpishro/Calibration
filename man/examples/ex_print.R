# build a simple calibration model and print the model estimates
cal1    <- calibrate(sim = Ds1, field = Df1,
                     Nmcmc = 5, nBurn = 0, thinning = 1,
                     kappaDist = "gaussian", kappaInit = 0.5, kappaP1 = 0.5, kappaP2 = 0.25)
print(cal1)
