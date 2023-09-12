# build a simple calibration model and print the model estimates
cal1    <- calibrate(sim = analytic11S, field = analytic11F,
                     nMCMC = 5, nBurn = 0, thinning = 1,
                     kappaDist = "normalTr", kappaP1 = 0.5, kappaP2 = 0.25)
print(cal1)
