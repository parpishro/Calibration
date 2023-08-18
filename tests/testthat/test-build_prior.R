test_that("build_prior produces correct functions", {
  expect_no_error(build_prior(dist = "gamma", p1 = 5, p2 = 2))
})
