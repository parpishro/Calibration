cal1    <- calibrate(sim = analytic11S, field = analytic11F,
                     nMCMC = 5, nBurn = 0, thinning = 1,
                     kappaDist = "normalTr", kappaP1 = 0.5, kappaP2 = 0.25)
summary(cal1)

cal2    <- calibrate(sim = Ds2, field = Df2, nMCMC = 20, nBurn = 0, thinning = 1)
summary(cal2)
