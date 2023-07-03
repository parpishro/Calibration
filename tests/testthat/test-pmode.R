test_that("multiplication works", {
  vec     <- c(1, 1, 1, 1, 1, 1, 1,
               2, 2, 2, 2,
               3,
               4, 4, 4, 4, 4, 4, 4, 4,
               5, 5,
               6, 6, 6,
               7, 7, 7, 7, 7, 7, 7, 7, 7,
               8, 8, 8,
               9, 9, 9)
  expect_equal(pmode(vec), 7)
})
