cal1    <- calibrate(sim = Ds1, field = Df1,
                     Nmcmc = 5, nBurn = 0, thinning = 1,
                     kappa = list(dist = "gaussian", init = 0.5, p1 = 0.5, p2 = 0.25))
summary(cal1)

cal2    <- calibrate(sim = Ds2, field = Df2, Nmcmc = 20, nBurn = 0, thinning = 1)
summary(cal2)
