test_that("prior_builder produces correct functions", {
  expect_no_error(prior_builder(prior = "gamma", p1 = 5, p2 = 2))
})
