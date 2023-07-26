test_that("plotting function works for simple calibration model", {
  cal   <- calibrate(sim = Ds1, field = Df1, nMCMC = 10, nBurn = 0, thinning = 1)
  expect_no_error(plot(cal))
  expect_no_error(plot(cal, parameter = "thetaS"))
  expect_no_error(plot(cal, parameter = "sigma2B", type = "trace"))
  expect_no_error(plot(cal, type = "fits"))
})

test_that("plotting function works for more complex calibration models", {
  cal   <- calibrate(sim = Ds2, field = Df2, nMCMC = 20, nBurn = 0, thinning = 1)
  expect_no_error(plot(cal))
  expect_no_error(plot(cal, parameter = "thetaS"))
  expect_no_error(plot(cal, parameter = "sigma2B", type = "trace"))
  expect_no_error(plot(cal, type = "fits"))
})
