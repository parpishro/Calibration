X     <- matrix(c(1, 3, 5,
                  2, 2, 6,
                  1, 4, 1), nrow = 3, byrow = TRUE)

Y     <- matrix(c(7, 3, 0,
                  2, 2, 4), nrow = 2, byrow = TRUE)

sc    <- c(1, 2, 3) # scale parameters of correlation structure
sm    <- c(2, 1, 2) # smoothness parameters of correlation structure

# correlation of a matrix with itself
correlation(X, theta = sc, alpha = sm)

# correlation between two matrices
correlation(X, Y, theta = sc, alpha = sm)

