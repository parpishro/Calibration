


## Load field data

swField  <- as.matrix(read.table("data-raw/spotweld/sw_field.txt", header = T)[-1])
usethis::use_data(swField)


## Load simulated data

swSim  <- as.matrix(na.omit(read.table("data-raw/spotweld/sw_sim.txt", header = T)[c(2, 4, 5, 3, 6)]))
usethis::use_data(swSim)


