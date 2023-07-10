test_that("Error in correlation of a matrix with itself!", {
  X     <- matrix(c(1, 3, 5,
                    2, 2, 6,
                    1, 4, 1), nrow = 3, byrow = TRUE)

  sc    <- c(1, 2, 3)
  sm    <- c(2, 1, 2)

  ExpXX <- matrix(c(1,        exp(-6),  exp(-50),
                    exp(-6),  1,        exp(-80),
                    exp(-50), exp(-80), 1), nrow = 3, byrow = TRUE)
  expect_equal(correlation(X, theta = sc, alpha = sm), ExpXX)
})

test_that("Error in correlation between two different matrices!", {
    X     <- matrix(c(1, 3, 5,
                    2, 2, 6,
                    1, 4, 1), nrow = 3, byrow = TRUE)

  Y     <- matrix(c(7, 3, 0,
                    2, 2, 4), nrow = 2, byrow = TRUE)

  sc    <- c(1, 2, 3)
  sm    <- c(2, 1, 2)

  ExpXY <- matrix(c(exp(-111), exp(-6),
                    exp(-135), exp(-12),
                    exp(-41),  exp(-32)), nrow = 3, byrow = TRUE)
  expect_equal(correlation(X, Y, theta = sc, alpha = sm), ExpXY)
})
