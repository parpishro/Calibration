test_that("MAP predict returns fbc object with correct elements  for 1 experimental
          + 1 calibration inputs", {
  cal   <- calibrate(sim = Ds1, field = Df1, nMCMC = 10, nBurn = 0, thinning = 1)
  preds <- predict(cal, newdata = matrix(c(2), nrow = 1), method = "MAP")
  expect_named(preds, expected = c("pred", "se"))
  types <- c(typeof(preds$pred), typeof(preds$se))
  expect_setequal(types, c("double", "double"))


})

test_that("MAP predict output elements have coorect lenght and contain correct elements
           for 1 experimental + 1 calibration inputs", {
  cal   <- calibrate(sim = Ds1, field = Df1, nMCMC = 10, nBurn = 0, thinning = 1)
  preds <- predict(cal, newdata = matrix(c(2), nrow = 1), method = "MAP")
  expect_equal(sum(!is.finite(preds$pred)), 0)
  expect_equal(sum(!is.finite(preds$se)), 0)
  expect_equal(length(preds$pred), 1)
  expect_equal(length(preds$se), 1)
})



test_that("MAP predict produces reasonable prediction and bounds for 1 experimental
          + 1 calibration inputs", {
  cal   <- calibrate(sim = Ds1, field = Df1, nMCMC = 10, nBurn = 0, thinning = 1)
  preds <- predict(cal, newdata = matrix(c(2), nrow = 1), method = "MAP")
  expect_true(preds$pred >= min(c(Ds1[, 1], Df1[, 1])))
  expect_true(preds$pred <= max(c(Ds1[, 1], Df1[, 1])))
  expect_true(preds$se > 0)
  expect_true(is.finite(preds$se))
})



test_that("Bayesian predict returns fbc object with correct elements  for 1 experimental
          + 1 calibration inputs", {
  cal   <- calibrate(sim = Ds1, field = Df1, nMCMC = 10, nBurn = 0, thinning = 1)
  preds <- predict(cal, newdata = matrix(c(2), nrow = 1), method = "Bayesian")
  expect_named(preds, expected = c("pred", "se"))
  types <- c(typeof(preds$pred), typeof(preds$se))
  expect_setequal(types, c("double", "double"))


})

test_that("Bayesian predict output elements have coorect lenght and contain correct
          elements  for 1 experimental + 1 calibration inputs", {
  cal   <- calibrate(sim = Ds1, field = Df1, nMCMC = 10, nBurn = 0, thinning = 1)
  preds <- predict(cal, newdata = matrix(c(2), nrow = 1), method = "Bayesian")
  expect_equal(sum(!is.finite(preds$pred)), 0)
  expect_equal(sum(!is.finite(preds$se)), 0)
  expect_equal(length(preds$pred), 1)
  expect_equal(length(preds$se), 1)
})



test_that("Bayesian predict produces reasonable prediction and bounds for 1 experimental
          + 1 calibration inputs", {
  cal   <- calibrate(sim = Ds1, field = Df1, nMCMC = 10, nBurn = 0, thinning = 1)
  preds <- predict(cal, newdata = matrix(c(2), nrow = 1), method = "Bayesian")
  expect_true(preds$pred >= min(c(Ds1[, 1], Df1[, 1])))
  expect_true(preds$pred <= max(c(Ds1[, 1], Df1[, 1])))
  expect_true(preds$se > 0)
  expect_true(is.finite(preds$se))
})


test_that("MAP predict returns fbc object with correct elements  for 2 experimental
          + 2 calibration inputs", {
  cal   <- calibrate(sim = Ds2, field = Df2, nMCMC = 10, nBurn = 0, thinning = 1)
  preds <- predict(cal, newdata = matrix(c(2, 2), nrow = 1), method = "MAP")
  expect_named(preds, expected = c("pred", "se"))
  types <- c(typeof(preds$pred), typeof(preds$se))
  expect_setequal(types, c("double", "double"))


})

test_that("MAP predict output elements have coorect lenght and contain correct elements
           for 2 experimental + 2 calibration inputs", {
  cal   <- calibrate(sim = Ds2, field = Df2, nMCMC = 10, nBurn = 0, thinning = 1)
  preds <- predict(cal, newdata = matrix(c(2, 2), nrow = 1), method = "MAP")
  expect_equal(sum(!is.finite(preds$pred)), 0)
  expect_equal(sum(!is.finite(preds$se)), 0)
  expect_equal(length(preds$pred), 1)
  expect_equal(length(preds$se), 1)
})



test_that("MAP predict produces reasonable prediction and bounds for 2 experimental
          + 2 calibration inputs", {
  cal   <- calibrate(sim = Ds2, field = Df2, nMCMC = 10, nBurn = 0, thinning = 1)
  preds <- predict(cal, newdata = matrix(c(2, 2), nrow = 1), method = "MAP")
  expect_true(preds$pred >= min(c(Ds2[, 1], Df2[, 1])))
  expect_true(preds$pred <= max(c(Ds2[, 1], Df2[, 1])))
  expect_true(preds$se > 0)
  expect_true(is.finite(preds$se))
})



test_that("Bayesian predict returns fbc object with correct elements  for 2 experimental
          + 2 calibration inputs", {
  cal   <- calibrate(sim = Ds2, field = Df2, nMCMC = 10, nBurn = 0, thinning = 1)
  preds <- predict(cal, newdata = matrix(c(2, 2), nrow = 1), method = "Bayesian")
  expect_named(preds, expected = c("pred", "se"))
  types <- c(typeof(preds$pred), typeof(preds$se))
  expect_setequal(types, c("double", "double"))


})

test_that("Bayesian predict output elements have coorect lenght and contain correct
          elements  for 2 experimental + 2 calibration inputs", {
  cal   <- calibrate(sim = Ds2, field = Df2, nMCMC = 10, nBurn = 0, thinning = 1)
  preds <- predict(cal, newdata = matrix(c(2, 2), nrow = 1), method = "Bayesian")
  expect_equal(sum(!is.finite(preds$pred)), 0)
  expect_equal(sum(!is.finite(preds$se)), 0)
  expect_equal(length(preds$pred), 1)
  expect_equal(length(preds$se), 1)
})



test_that("Bayesian predict produces reasonable prediction and bounds for 2 experimental
          + 2 calibration inputs", {
  cal   <- calibrate(sim = Ds2, field = Df2, nMCMC = 10, nBurn = 0, thinning = 1)
  preds <- predict(cal, newdata = matrix(c(2, 2), nrow = 1), method = "Bayesian")
  expect_true(preds$pred >= min(c(Ds2[, 1], Df2[, 1])))
  expect_true(preds$pred <= max(c(Ds2[, 1], Df2[, 1])))
  expect_true(preds$se > 0)
  expect_true(is.finite(preds$se))
})

test_that("Bayesian predict output elements have coorect lenght and contain correct
          elements  for 2 experimental + 2 calibration inputs (multiple inputs)", {
  cal   <- calibrate(sim = Ds2, field = Df2, nMCMC = 10, nBurn = 0, thinning = 1)
  preds <- predict(cal, newdata = matrix(c(2, 2, 3, 4), nrow = 2), method = "Bayesian")
  expect_equal(sum(!is.finite(preds$pred)), 0)
  expect_equal(sum(!is.finite(preds$se)), 0)
  expect_equal(length(preds$pred), 2)
  expect_equal(length(preds$se), 2)
})


