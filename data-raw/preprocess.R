library(lhs)


## Load field data

toyField  <- as.matrix(read.csv("data-raw/ball.csv"))
usethis::use_data(toyField)


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
sim <- function(XT) sqrt(2*XT[, 1]/XT[, 2])

# run simulation (based on field ranges)

toySim   <- matrix(nrow = 100, ncol = 3)
toySim[, 1:2]   <- design_lhs(m      = 100,    # number of simulator runs
                              inputs = 2,      # p+q = 2
                              ranges = list(h = range(toyField[, 'height']),
                                            g = c(6, 14)))
toySim[, 3] <- sim(toySim[, 1:2])  + rnorm(100, mean = 0, sd = 0.01) # simulated output (time)

# Save the cleaned data in the required R package location
usethis::use_data(toySim)


# parameters

toyParams <- list(theta0 = 9.8, omega0 = 0.1, alpha0 = 5, sigma20 = 1)
