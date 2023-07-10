create_data <- function(n, m, p, q, sdE, fnM, fnS, seed = 666) {
  set.seed(seed)
  Xf  <- lhs::maximinLHS(n, p) * 10 - 5
  Xs  <- lhs::maximinLHS(m, p + q)
  scl <- c(rep(10, p), rep(4, q))
  Xs  <- t(apply(Xs, 1, function(x) x* scl - (scl/2)))
  err <- rnorm(n, 0, sdE)
  yf  <- fnM(Xf) + err
  ys  <- fnS(Xs)
  Df  <- cbind(matrix(yf, ncol = 1), Xf)
  Ds  <- cbind(matrix(ys, ncol = 1), Xs)
  return(list(Ds = Ds, Df = Df))
}

fnM1 <- function(X) 0.1*(X[,1]^2) - X[1,] + 0.4
fnS1 <- function(X) X[,2]*(X[,1]^2) - 1.3*X[, 1]
D1   <- create_data(n = 5, m = 15, p = 1, q = 1, sdE = 0.8, fnM = fnM1, fnS = fnS1)
Df1  <- D1$Df
Ds1  <- D1$Ds

fnM2 <- function(X) 0.1*(X[,1]^2) + X[, 2]
fnS2 <- function(X) X[, 3] * (X[, 1] ^ 2) + X[, 4] * X[, 2] + 1.1
D2   <- create_data(n = 50, m = 100, p = 2, q = 2, sdE = 0.24, fnM = fnM2, fnS = fnS2)
Df2  <- D2$Df
Ds2  <- D2$Ds

usethis::use_data(Df1)
usethis::use_data(Ds1)
usethis::use_data(Df2)
usethis::use_data(Ds2)
