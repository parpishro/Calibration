cal     <- calibrate(sim = Ds, field = Df,
                     Nmcmc = 5L, nBurn = 0L, thinning = 1L,
                     kappa_dist = "gaussian", init = 0.5, p1 = 0.5, p2 = 0.25)
summary(cal)
