#' compute the correlation matrix of between the rows of two matrices.
#'
#' If only one matrix is provided, `correlation` would compute the correlation between
#' the rows of the given matrix. In this formulation, the output correlation matriz will be symmetric.
#' Power exponential correlation structure is assumed. Each entry R[i, j] in the returned
#' correlation matrix is correlation between row i of X and row j of Y (X if Y is null) using
#' smaoothness and scale parameters. Distance between rows, and smoothness and scale
#  parameters are computed element-by-element using Euclidean metric.
#'
#' @param X      A matrix with at least two rows
#' @param Y      A matrix with same number of columns as X
#' @param theta  A vector of scale parameters ($\theta \in (0, \inf)$) with the same size as number of columns as X
#' @param alpha  A vector of smoothness parameters ($\alpha \in [1, 2]$) with the same size as number of columns as X
#'
#' @return a correlation matrix between the rows of given matrices
correlation <- function (X, Y = NULL, theta, alpha) {
  nx    <- nrow(X)
  if (is.null(Y)) {
    R <- matrix(0, nrow = nx, ncol = nx)
    for (i in 1:nx-1) {
      for (j in (i+1):nx) {
        R[i, j] <- exp(-sum(theta * (abs(X[i, ] - X[j, ]) ^ alpha)))
      }
    }
    R <- R + t(R) + diag(1, nrow = nx, ncol = nx)
  } else {
    ny <- nrow(Y)
    R  <- matrix(0, nrow = nx, ncol = ny)
    for (i in 1:nx) {
      for (j in 1:ny) {
        R[i, j] <- exp(-sum(theta* (abs(X[i, ] - Y[j, ]) ^ alpha)))
      }
    }
  }
  return(R)
}

