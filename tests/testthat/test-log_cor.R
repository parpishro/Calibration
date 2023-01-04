X       <- matrix(c(1, 3, 5,
                    2, 2, 6,
                    1, 4, 1), nrow = 3, byrow = TRUE)

Y       <- matrix(c(7, 3, 0,
                    2, 2, 4), nrow = 2, byrow = TRUE)

scale   <- c(1, 2, 3)
smooth  <- c(2, 1, 2)

ExpXX   <- matrix(c(0,  -12,  -146,
                  -12,    0,    -230,
                 -146,   -230, 0), nrow = 3, byrow = TRUE)

ExpXY   <- matrix(c(-261, -12,
                    -351, -36,
                    -47,  -86), nrow = 3, byrow = TRUE)



test_that("distance with self works!", {
  expect_equal(log_cor(X, scale = scale, smooth = smooth), ExpXX)
})

test_that("distance of two matrices work!", {
  expect_equal(log_cor(X, Y, scale = scale, smooth = smooth), ExpXY)
})

