library(lhs)


## Load field data

ballField           <- as.matrix(read.csv("data-raw/ball/ball.csv")[c(2, 1)])
colnames(ballField) <- c("t", "h")
usethis::use_data(ballField)


## Create simulated data

# design input matrix
design_lhs <- function(m, inputs, ranges) {
  D <- maximinLHS(m, inputs)
  for (i in 1:inputs) {
    D[ ,i] <- ranges[[i]][1] + (D[ ,i] * diff(ranges[[i]]))
  }
  return(D)
}

# simulation process
sim <- function(X) sqrt(2*X[, 1]/X[, 2])

# run simulation (based on field ranges)

ballSim   <- matrix(nrow = 100, ncol = 3)
ballSim[, 2:3]   <- design_lhs(m      = 100,    # number of simulator runs
                               inputs = 2,      # p+q = 2
                               ranges = list(h = range(ballField[, 'h']),
                                             g = c(6, 14)))
ballSim[, 1] <- sim(ballSim[, 2:3])             # simulated output (time)
ballSim      <- round(ballSim, 3)
colnames(ballSim) <- c("t", "h", "g")
# Save the cleaned data in the required R package location

usethis::use_data(ballSim)

