# compute the correlation matrix of a matrix or two matrices

# REQUIRES: - X and Y must be matrix and scale and smoothness must be vector
#           - X must have at least two rows
#           - If Y is provided, it must have same number of columns as X
#           - scale and smoothness must have same number of elements as columns
#             in X (and Y)
# EFFECTS: computes the correlation matrix between two given matrices X and Y.
#             if Y is null, it computes the correlation matrix of X. Power
#             exponential correlation structure is assumed. Each entry R[i, j]
#             in the returned correlation matrix is correlation between row i
#             of X and row j of Y (X if Y is null) using smaoothness and scale
#             parameters. Distance between rows, and smoothness and scale
#             parameters are computed element-by-element, When Y is null, a
#             branch makes the computation easier (R is symmetric).

correlation <- function (X, Y = NULL, lambda, gamma) {
  nx    <- nrow(X)
  rho   <- lambda
  alpha <- 1 + 1/(1+exp(-gamma))
  if (is.null(Y)) {
    R <- matrix(0, nrow = nx, ncol = nx)
    for (i in 1:nx-1) {
      for (j in (i+1):nx) {
        R[i, j] <- prod(rho) * exp(-sum(abs(X[i, ] - X[j, ]) ^ alpha))
      }
    }
    R <- R + t(R) + diag(1, nrow = nx, ncol = nx)
  } else {
    ny <- nrow(Y)
    R  <- matrix(0, nrow = nx, ncol = ny)
    for (i in 1:nx) {
      for (j in 1:ny) {
        R[i, j] <- prod(rho) * exp(-sum(abs(X[i, ] - Y[j, ]) ^ alpha))
      }
    }
  }
  return(R)
}

