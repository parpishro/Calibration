test_that("print produces output in console", {
  cal   <- calibrate(sim = analytic11S, field = analytic11F,
                     nMCMC = 10, nBurn = 0, thinning = 1)
  expect_invisible(print(cal))
  expect_output(print(cal))
})
