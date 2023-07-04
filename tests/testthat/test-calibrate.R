test_that("calibrate returns fbc object with correct elements", {
  cal <- calibrate(sim = Ds1, field = Df1, Nmcmc = 10, nBurn = 0, thinning = 1)
  expect_s3_class(cal, "fbc")
  expect_named(cal, expected = c("Phi", "estimates", "logPost", "priors", "acceptance", "vars", "cache"))
  types <- c(typeof(cal$Phi), typeof(cal$estimates), typeof(cal$logPost), typeof(cal$priors),
             typeof(cal$acceptance), typeof(cal$vars), typeof(cal$cache))
  expect_setequal(types, c("list", "list", "double", "list", "double", "character", "environment"))


})

test_that("calibrate output elements have coorect lenght and contain correct elements", {
  cal <- calibrate(sim = Ds1, field = Df1, Nmcmc = 10, nBurn = 0, thinning = 1)
  expect_equal(sum(!is.finite(as.matrix(cal$Phi))), 0)
  expect_equal(sum(!is.finite(as.matrix(cal$estimates))), 0)
  expect_equal(sum(!is.finite(cal$logPost)), 0)
  expect_equal(sum(!is.finite(cal$acceptance)), 0)
  expect_named(cal$estimates, expected = c("mean", "median", "mode", "lwr50", "upr50", "lwr80", "upr80", "sd"))
  labels <- c(paste0("kappa", 1:cal$cache$q), paste0("thetaS", 1:(cal$cache$p + cal$cache$q)),
              paste0("alphaS", 1:(cal$cache$p + cal$cache$q)),paste0("thetaB", 1:cal$cache$p),
              paste0("alphaB", 1:cal$cache$p), "sigma2S", "sigma2B", "sigma2E", "muB")
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
  c   <- cal$cache
  for (i in c$ikappa) {
    expect_equal(sum(cal$estimates[i, -8] < min(Ds1[,1 + c$p+i])), 0)
    expect_equal(sum(cal$estimates[i, -8] > max(Ds1[,1 + c$p+i])), 0)
  }
  for (i in c(c$ithetaS, c$ithetaB, c$isigma2S, c$isigma2B, c$isigma2E))
    expect_equal(sum(cal$estimates[i, -8] < 0), 0)
  for (i in c(c$ialphaS, c$ialphaB)) {
    expect_equal(sum(cal$estimates[i, -8] < 1), 0)
    expect_equal(sum(cal$estimates[i, -8] > 2), 0)
  }
  expect_equal(sum(cal$estimates[c$imuB, -8] < -1), 0)
  expect_equal(sum(cal$estimates[c$imuB, -8] > 1), 0)

})


test_that("adaptive proposal results in proper mixing!", {
  cal <- calibrate(sim = Ds1, field = Df1, Nmcmc = 10, nBurn = 0, thinning = 1)
  c   <- cal$cache
  expect_equal(sum(!is.finite(cal$estimates[, 8])), 0)
  expect_equal(sum(cal$estimates[, 8] <= 0), 0)
  expect_equal(sum(cal$acceptance < 0.1), 0)
  expect_equal(sum(cal$acceptance > 0.9), 0)
})


test_that("calibrate returns fbc object with correct elements", {
  cal <- calibrate(sim = Ds2, field = Df2, Nmcmc = 10, nBurn = 0, thinning = 1)
  expect_s3_class(cal, "fbc")
  expect_named(cal, expected = c("Phi", "estimates", "logPost", "priors", "acceptance", "vars", "cache"))
  types <- c(typeof(cal$Phi), typeof(cal$estimates), typeof(cal$logPost), typeof(cal$priors),
             typeof(cal$acceptance), typeof(cal$vars), typeof(cal$cache))
  expect_setequal(types, c("list", "list", "double", "list", "double", "character", "environment"))


})

test_that("calibrate output elements have coorect lenght and contain correct elements", {
  cal <- calibrate(sim = Ds2, field = Df2, Nmcmc = 10, nBurn = 0, thinning = 1)
  expect_equal(sum(!is.finite(as.matrix(cal$Phi))), 0)
  expect_equal(sum(!is.finite(as.matrix(cal$estimates))), 0)
  expect_equal(sum(!is.finite(cal$logPost)), 0)
  expect_equal(sum(!is.finite(cal$acceptance)), 0)
  expect_named(cal$estimates, expected = c("mean", "median", "mode", "lwr50", "upr50", "lwr80", "upr80", "sd"))
  labels <- c(paste0("kappa", 1:cal$cache$q), paste0("thetaS", 1:(cal$cache$p + cal$cache$q)),
              paste0("alphaS", 1:(cal$cache$p + cal$cache$q)),paste0("thetaB", 1:cal$cache$p),
              paste0("alphaB", 1:cal$cache$p), "sigma2S", "sigma2B", "sigma2E", "muB")
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
  c   <- cal$cache
  for (i in c$ikappa) {
    expect_equal(sum(cal$estimates[i, -8] < min(Ds2[,1 + c$p+i])), 0)
    expect_equal(sum(cal$estimates[i, -8] > max(Ds2[,1 + c$p+i])), 0)
  }
  for (i in c(c$ithetaS, c$ithetaB, c$isigma2S, c$isigma2B, c$isigma2E))
    expect_equal(sum(cal$estimates[i, -8] < 0), 0)
  for (i in c(c$ialphaS, c$ialphaB)) {
    expect_equal(sum(cal$estimates[i, -8] < 1), 0)
    expect_equal(sum(cal$estimates[i, -8] > 2), 0)
  }
  expect_equal(sum(cal$estimates[c$imuB, -8] < -1), 0)
  expect_equal(sum(cal$estimates[c$imuB, -8] > 1), 0)

})


test_that("adaptive proposal results in proper mixing!", {
  cal <- calibrate(sim = Ds2, field = Df2, Nmcmc = 10, nBurn = 0, thinning = 1)
  c   <- cal$cache
  expect_equal(sum(!is.finite(cal$estimates[, 8])), 0)
  expect_equal(sum(cal$estimates[, 8] <= 0), 0)
  expect_equal(sum(cal$acceptance < 0.1), 0)
  expect_equal(sum(cal$acceptance > 0.9), 0)
})
