X <- matrix(c(1, 3, 5,
              2, 2, 6,
              1, 4, 1), nrow = 3, byrow = TRUE)

Y <- matrix(c(7, 3, 0,
              2, 2, 4), nrow = 2, byrow = TRUE)

ExpXX <- matrix(c(0, 0, 0,
                  1, 1, 1,
                  0, 1, 4,
                  1, 1, 1,
                  0, 0, 0,
                  1, 2, 5,
                  0, 1, 4,
                  1, 2, 5,
                  0, 0, 0), nrow = 9, byrow = TRUE)
ExpXY <- matrix(c(6, 0, 5,
                  1, 1, 1,
                  5, 1, 6,
                  0, 0, 2,
                  6, 1, 1,
                  1, 2, 3), nrow = 6, byrow = TRUE)



test_that("distance with self works!", {
  expect_equal(distance(X, X), ExpXX)
})

test_that("distance of two matrices work!", {
  expect_equal(distance(X, Y), ExpXY)
})



