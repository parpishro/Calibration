#' Correlation Between Rows of Two Matrices.
#'
#' If only one matrix is provided, `correlation` would compute the correlation between
#' the rows of the given matrix. Each entry \eqn{R[i, j]} in the returned correlation matrix is
#' correlation between row i of X and row j of Y (X if Y is null) using given smoothness
#' and scale parameters. To compute correlation, separable power exponential correlation
#' structure is assumed and using Euclidean notion of distance between rows.
#'
#' @param X      matrix with at least two rows
#' @param Y      matrix with same number of columns as X
#' @param theta  vector of scale parameters (\eqn{\theta \in (0, \inf)}) with the same
#'               size as number of columns as X
#' @param alpha  vector of smoothness parameters (\eqn{\alpha \in [1, 2]}) with the same
#'               size as number of columns as X
#'
#' @returns      a correlation matrix between the rows of given matrices
#' @example      man/examples/ex_correlation.R
#' @export
correlation <- function(X, Y = NULL, theta, alpha) {
  nx <- nrow(X)
  if (is.null(Y)) {
    R  <- matrix(0, nrow = nx, ncol = nx)
    for (i in 1:(nx - 1)) {
      for (j in (i + 1):nx) {
        R[i, j] <- -sum(theta*(abs(X[i, ] - X[j, ])^alpha))
      }
    }
    R <- R + t(R) + diag(0, nrow = nx, ncol = nx)
  } else {
    ny <- nrow(Y)
    R  <- matrix(0, nrow = nx, ncol = ny)
    for (i in 1:nx) {
      for (j in 1:ny) {
        R[i, j] <- -sum(theta*(abs(X[i, ] - Y[j, ])^alpha))
      }
    }
  }
  return(exp(R))
}





