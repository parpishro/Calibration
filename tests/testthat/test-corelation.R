X       <- matrix(c(1, 3, 5,
                    2, 2, 6,
                    1, 4, 1), nrow = 3, byrow = TRUE)

Y       <- matrix(c(7, 3, 0,
                    2, 2, 4), nrow = 2, byrow = TRUE)

scale   <- c(1, 2, 3)
smooth  <- c(2, 1, 2)

DXY <- matrix(c(6, 0, 5,
                1, 1, 1,
                5, 1, 6,
                0, 0, 2,
                6, 1, 1,
                1, 2, 3), nrow = 6, byrow = TRUE)

ExpXX   <- matrix(c(1,         exp(-6),   exp(-50),
                    exp(-6),   1,         exp(-80),
                    exp(-50),  exp(-80),  1), nrow = 3, byrow = TRUE)

ExpXY   <- matrix(c(exp(-111), exp(-6),
                    exp(-135), exp(-12),
                    exp(-41),  exp(-32)), nrow = 3, byrow = TRUE)



test_that("log_cor with self works!", {
  expect_equal(corelation(X, scale = scale, smooth = smooth), ExpXX)
})

test_that("log_cor of two matrices work!", {
  expect_equal(corelation(X, Y, scale = scale, smooth = smooth), ExpXY)
})

