cal1    <- calibrate(sim = Ds1, field = Df1,
                     Nmcmc = 5L, nBurn = 0L, thinning = 1L,
                     kappa_dist = "gaussian", init = 0.5, p1 = 0.5, p2 = 0.25)
summary(cal1)

cal2    <- calibrate(sim = Ds2, field = Df2, Nmcmc = 20L, nBurn = 0L, thinning = 1L)
summary(cal2)
