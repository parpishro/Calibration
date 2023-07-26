test_that("summary works as expected for simple calibration", {
  cal   <- calibrate(sim = Ds1, field = Df1, nMCMC = 10, nBurn = 0, thinning = 1)
  expect_no_error(summary(cal))
})


test_that("summary works as expected for complex calibration", {
  cal   <- calibrate(sim = Ds2, field = Df2, nMCMC = 20, nBurn = 0, thinning = 1)
  expect_no_error(summary(cal))
})
