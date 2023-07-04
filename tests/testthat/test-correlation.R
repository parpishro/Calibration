test_that("Error in correlation of a matrix with itself!", {
  dat <- create_testData()
  expect_equal(correlation(dat$X, theta = dat$sc, alpha = dat$sm), dat$ExpXX)
})

test_that("Error in correlation between two different matrices!", {
  dat <- create_testData()
  expect_equal(correlation(dat$X, dat$Y, theta = dat$sc, alpha = dat$sm), dat$ExpXY)
})
