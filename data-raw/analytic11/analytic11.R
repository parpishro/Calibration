library(lhs)

xf          <- c(-5, -2.5, 0, 2.5, 5)
yf          <- c(7.079, 4.538, -1.539, -2.096, -0.824) # yf(xf) = 0.1x^2 - x + 0.4 + N(0, 0.8^2)
analytic11F <- as.matrix(cbind(yf, xf))

Xs          <- (maximinLHS(15, 2) - 0.5) * 10
ys          <- Xs[, 2] - 1.3*Xs[, 1]
analytic11S <- as.matrix(cbind(ys, Xs))



usethis::use_data(analytic11F)
usethis::use_data(analytic11S)



set.seed(666)
n         <- 30
p         <- 2
q         <- 2
kappa     <- c(0.25,0.4)
sdE       <- 0.1

Xf        <- 3*lhs::maximinLHS(n, p)
Xk        <- cbind(Xf, matrix(rep(kappa, n), ncol = q, nrow = n, byrow = T))
Eps       <- rnorm(n, 0, sdE)

fnS       <- function(X) -((X[, 3]/6)*(X[, 2]^2)) + (sqrt(X[, 4])*sqrt(X[,2])) + (0.1*X[, 1])
fnB       <- function(X) 0.1*(X[,1]^2) - (0.2*(X[,2]*sqrt(X[,2])))

yf        <- fnS(Xk) + fnB(Xf) + Eps

DfFull    <- as.matrix(data.frame(yf   = yf,
                                  x1   = Xf[,1],
                                  x2   = Xf[,2],
                                  sim  = fnS(Xk),
                                  bias = fnB(Xf),
                                  eps  = Eps))

Df2       <- DfFull[, 1:3]

m         <- 60
Xs        <- cbind(3*lhs::maximinLHS(m, p), lhs::maximinLHS(m, q))
ys        <- fnS(Xs)
Ds2       <- as.matrix(data.frame(ys = ys,
                                  x1 = Xs[,1],
                                  x2 = Xs[,2],
                                  k1 = Xs[,3],
                                  k2 = Xs[,4]))

usethis::use_data(Df2)
usethis::use_data(Ds2)
