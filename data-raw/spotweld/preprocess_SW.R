


## Load field data

swField  <- as.matrix(read.table("data-raw/spotweld/sw_field.txt", header = T)[c(5,2,3,4)])
usethis::use_data(swField)


## Load simulated data

swSim  <- as.matrix(na.omit(read.table("data-raw/spotweld/sw_sim.txt", header = T)[c(6, 2, 4, 5, 3)]))
usethis::use_data(swSim)


