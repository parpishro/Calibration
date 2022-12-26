# compute the correlation matrix of a matrix or two matrices

# REQUIRES: X and Y must be matrix and cale and smoothness must be vector. If Y
#               is provided, it must have same number of columns as X. length
#               scale and smoothness vectors must be same as number of columns
#               in X (and Y)
# EFFECTS:

correlation <- function (X, Y = NULL, scale, smoothness) {
  if (is.null(Y)) Y <- X
  D <- distance(X, Y)
  nx <- nrow(X)
  ny <- nrow(Y)
  R <- matrix(nrow = nx, ncol = ny)
  for (i in 1:nx) {
    for (j in 1:ny) {
      R[i, j] <- prod(exp(-(D[((i - 1) * ny) + j, ] * scale) ^ smoothness))
    }
  }
  return(R)
}
