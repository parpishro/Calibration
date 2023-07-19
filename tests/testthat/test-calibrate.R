test_that("calibrate returns fbc object with correct elements", {
  cal <- calibrate(sim = Ds1, field = Df1, Nmcmc = 10, nBurn = 0, thinning = 1)
  expect_s3_class(cal, "fbc")
  expect_named(cal, expected = c("Phi", "estimates", "logPost", "priors", "acceptance", "vars",
                                 "data", "scale", "indices", "priorFns", "proposalSD"))
  types <- c(typeof(cal$Phi), typeof(cal$estimates), typeof(cal$logPost), typeof(cal$priors),
             typeof(cal$acceptance), typeof(cal$vars), typeof(cal$data), typeof(cal$scale),
             typeof(cal$indices), typeof(cal$priorFns), typeof(cal$proposalSD))
  expect_setequal(types, c("list", "list", "double", "list", "double", "character", "list", "list",
                           "list", "list", "double"))


})

test_that("calibrate output elements have coorect lenght and contain correct elements", {
  cal <- calibrate(sim = Ds1, field = Df1, Nmcmc = 10, nBurn = 0, thinning = 1)
  expect_equal(sum(!is.finite(as.matrix(cal$Phi))), 0)
  expect_equal(sum(!is.finite(as.matrix(cal$estimates))), 0)
  expect_equal(sum(!is.finite(cal$logPost)), 0)
  expect_equal(sum(!is.finite(cal$acceptance)), 0)
  expect_named(cal$estimates, expected = c("mean", "median", "mode", "lwr50", "upr50", "lwr80", "upr80", "sd"))
  labels <- c(paste0("kappa", 1:cal$indices$q), paste0("thetaS", 1:(cal$indices$p + cal$indices$q)),
              paste0("alphaS", 1:(cal$indices$p + cal$indices$q)),paste0("thetaB", 1:cal$indices$p),
              paste0("alphaB", 1:cal$indices$p), "sigma2S", "sigma2B", "sigma2E", "muB")
  expect_named(cal$Phi, expected = labels)
  expect_setequal(cal$vars, labels)
  expect_equal(nrow(cal$Phi), 10)
  expect_equal(nrow(cal$estimates), length(labels))
  expect_equal(length(cal$logPost), 10)
  expect_equal(length(cal$acceptance), length(labels))
  expect_equal(length(cal$priors), 9)

})



test_that("posterior parameter distribution produces reasonable summary statistics!", {
  cal <- calibrate(sim = Ds1, field = Df1, Nmcmc = 10, nBurn = 0, thinning = 1)
  for (i in cal$indices$ikappa) {
    expect_equal(sum(cal$estimates[i, -8] < min(Ds1[,1 + cal$indices$p + i])), 0)
    expect_equal(sum(cal$estimates[i, -8] > max(Ds1[,1 + cal$indices$p + i])), 0)
  }
  for (i in c(cal$indices$ithetaS, cal$indices$ithetaB, cal$indices$isigma2S,
              cal$indices$isigma2B, cal$indices$isigma2E))
    expect_equal(sum(cal$estimates[i, -8] < 0), 0)
  for (i in c(cal$indices$ialphaS, cal$indices$ialphaB)) {
    expect_equal(sum(cal$estimates[i, -8] < 1), 0)
    expect_equal(sum(cal$estimates[i, -8] > 2), 0)
  }
  expect_equal(sum(cal$estimates[cal$indices$imuB, -8] < -1), 0)
  expect_equal(sum(cal$estimates[cal$indices$imuB, -8] > 1), 0)

})


test_that("adaptive proposal results in proper mixing!", {
  cal <- calibrate(sim = Ds1, field = Df1, Nmcmc = 10, nBurn = 0, thinning = 1)
  expect_equal(sum(!is.finite(cal$estimates[, 8])), 0)
  expect_equal(sum(cal$estimates[, 8] <= 0), 0)
  expect_equal(sum(cal$acceptance == 0), 0)
})


test_that("fixed parameter works!", {
  pr  <- set_hyperPriors(alphaSDist = "fixed", alphaSInit = 2)
  cal <- calibrate(sim = Ds1, field = Df1, Nmcmc = 10, nBurn = 0, thinning = 1, hypers = pr)
  expect_equal(sum(!is.finite(cal$estimates[, 8])), 0)
  expect_equal(sum(cal$estimates[, 8] < 0), 0)
  expect_equal(sum(cal$acceptance[c(1:3, 6:11)] == 0), 0)
})


test_that("different priors for calibration works!", {
  cal <- calibrate(sim = Ds2, field = Df2, Nmcmc = 15, nBurn = 0, thinning = 1,
                   kappaDist = c("beta", "uniform"),
                                kappaInit = c(0.5, 0.3),
                                kappaP1   = c(1.1, 0),
                                kappaP2   = c(1.1, 1))
  expect_equal(sum(!is.finite(cal$estimates[, 8])), 0)
  expect_equal(sum(cal$estimates[, 8] <= 0), 0)
})


test_that("progressBar works!", {
  cal <- calibrate(sim = Ds1, field = Df1, Nmcmc = 200, nBurn = 0, thinning = 1, showProgress = TRUE)
  expect_equal(sum(!is.finite(cal$estimates[, 8])), 0)
  expect_equal(sum(cal$estimates[, 8] <= 0), 0)
  expect_equal(sum(cal$acceptance == 0), 0)
})


test_that("calibrate returns fbc object with correct elements", {
  cal <- calibrate(sim = Ds2, field = Df2, Nmcmc = 10, nBurn = 0, thinning = 1)
  expect_s3_class(cal, "fbc")
  expect_named(cal, expected = c("Phi", "estimates", "logPost", "priors", "acceptance", "vars",
                                  "data", "scale", "indices", "priorFns", "proposalSD"))
  types <- c(typeof(cal$Phi), typeof(cal$estimates), typeof(cal$logPost), typeof(cal$priors),
             typeof(cal$acceptance), typeof(cal$vars), typeof(cal$data), typeof(cal$scale),
             typeof(cal$indices), typeof(cal$priorFns), typeof(cal$proposalSD))
  expect_setequal(types, c("list", "list", "double", "list", "double", "character", "list", "list",
                           "list", "list", "double"))


})

test_that("calibrate output elements have coorect lenght and contain correct elements", {
  cal <- calibrate(sim = Ds2, field = Df2, Nmcmc = 10, nBurn = 0, thinning = 1)
  expect_equal(sum(!is.finite(as.matrix(cal$Phi))), 0)
  expect_equal(sum(!is.finite(as.matrix(cal$estimates))), 0)
  expect_equal(sum(!is.finite(cal$logPost)), 0)
  expect_equal(sum(!is.finite(cal$acceptance)), 0)
  expect_named(cal$estimates, expected = c("mean", "median", "mode", "lwr50", "upr50", "lwr80", "upr80", "sd"))
  labels <- c(paste0("kappa", 1:cal$indices$q), paste0("thetaS", 1:(cal$indices$p + cal$indices$q)),
              paste0("alphaS", 1:(cal$indices$p + cal$indices$q)),paste0("thetaB", 1:cal$indices$p),
              paste0("alphaB", 1:cal$indices$p), "sigma2S", "sigma2B", "sigma2E", "muB")
  expect_named(cal$Phi, expected = labels)
  expect_setequal(cal$vars, labels)
  expect_equal(nrow(cal$Phi), 10)
  expect_equal(nrow(cal$estimates), length(labels))
  expect_equal(length(cal$logPost), 10)
  expect_equal(length(cal$acceptance), length(labels))
  expect_equal(length(cal$priors), 9)

})



test_that("posterior parameter distribution produces reasonable summary statistics!", {
  cal <- calibrate(sim = Ds2, field = Df2, Nmcmc = 10, nBurn = 0, thinning = 1)
  for (i in cal$indices$ikappa) {
    expect_equal(sum(cal$estimates[i, -8] < min(Ds2[,1 + cal$indices$p + i])), 0)
    expect_equal(sum(cal$estimates[i, -8] > max(Ds2[,1 + cal$indices$p + i])), 0)
  }
  for (i in c(cal$indices$ithetaS, cal$indices$ithetaB, cal$indices$isigma2S, cal$indices$isigma2B,
              cal$indices$isigma2E))
    expect_equal(sum(cal$estimates[i, -8] < 0), 0)
  for (i in c(cal$indices$ialphaS, cal$indices$ialphaB)) {
    expect_equal(sum(cal$estimates[i, -8] < 1), 0)
    expect_equal(sum(cal$estimates[i, -8] > 2), 0)
  }
  expect_equal(sum(cal$estimates[cal$indices$imuB, -8] < -1), 0)
  expect_equal(sum(cal$estimates[cal$indices$imuB, -8] > 1), 0)

})


test_that("adaptive proposal results in proper mixing!", {
  cal <- calibrate(sim = Ds2, field = Df2, Nmcmc = 10, nBurn = 0, thinning = 1)
  expect_equal(sum(!is.finite(cal$estimates[, 8])), 0)
  expect_equal(sum(cal$estimates[, 8] <= 0), 0)
  expect_equal(sum(cal$acceptance == 0), 0)
})
