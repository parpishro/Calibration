library(lhs)

set.seed(666)
n   <- 5
m   <- 15
p   <- 1
q   <- 1
sdE <- 0.8
fnM <- function(x) 0.1*(x^2) - x + 0.4
err <- rnorm(n, 0, sdE) # error mean = 0, error sd = 0.8, n = 5
xf  <- seq(-5, 5, 2.5)
yf  <- fnM(xf) + err
Df  <- cbind(matrix(yf, ncol = 1), matrix(xf, ncol=1))



fnS <- function(x, k) k - 1.3*x
Xs  <- maximinLHS(m, p + q) * c(-5, 5) # m = 15, p + q = 2
ys  <- fnS(Xs[,1:p], Xs[,(p + 1):(p + q)])

Ds <- cbind(matrix(ys, ncol = 1), Xs)
Ds <- Ds[order(Ds[, 2]), ]

#ylimit <- c(min(c(yf, ys)), max(c(yf, ys)))
#plot(Df[, 2], Df[, 1], type = "l", ylim = ylimit)
#lines(seq(-5, 5, 0.01), fnM(seq(-5, 5, 0.01)), col = "blue")
#lines(Ds[, 2], Ds[, 1], col = "red")

usethis::use_data(Df)
usethis::use_data(Ds)
