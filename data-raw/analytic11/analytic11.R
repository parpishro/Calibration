library(lhs)

xf          <- c(-5, -2.5, 0, 2.5, 5)
yf          <- c(7.079, 4.538, -1.539, -2.096, -0.824) # yf(xf) = 0.1x^2 - x + 0.4 + N(0, 0.8^2)
analytic11F <- cbind(yf, xf)

Xs          <- (maximinLHS(15, 2) - 0.5) * 10
ys          <- Xs[, 2] - 1.3*Xs[, 1]
analytic11S <- cbind(ys, Xs)

usethis::use_data(analytic11F)
usethis::use_data(analytic11S)
