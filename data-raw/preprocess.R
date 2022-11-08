library(lhs)


## Load field data

field   <- read.csv("data-raw/ball.csv")
toy1_z  <- field$time
toy1_Xz <- matrix(field$height, ncol = 1)
# Save the cleaned data in the required R package location
usethis::use_data(toy1_z)
usethis::use_data(toy1_Xz)

## Create simulated data

# design input matrix
design_lhs <- function(m, inputs, ranges) {
  lhs_std <- maximinLHS(m, inputs)
  for (i in 1:inputs) {
    lhs_std[ ,i] <- ranges[[i]][1] + (lhs_std[ ,i] * diff(ranges[[i]]))
  }
  return(lhs_std)
}

# simulation process
sim <- function(XT) sqrt(2*XT[, 1]/XT[, 2])

# run simulation (based on field ranges)

toy1_XTy <- design_lhs(m      = 100,    # number of simulator runs
                       inputs = 2,     # d+q = 2
                       ranges = list(h = range(toy1_Xz),
                                     g = c(6, 14))) #expert opinion
toy1_y <- sim(toy1_XTy) # simulated output (time)

# Save the cleaned data in the required R package location
usethis::use_data(toy1_XTy)
usethis::use_data(toy1_y)
