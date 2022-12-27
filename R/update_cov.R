update_cov <- function() {
  if (changed == "calibration") {

  } else if (changed == "omega sim") {

  } else if (changed == "alpha sim") {

  } else if (changed == "omega bias") {

  } else if (changed == "alpha bias") {

  } else if (changed == "variance sim") {

  } else if (changed == "variance bias") {

  } else if (changed == "variance sim") {

  } else if (is.null(changed)) {

  } else stop("invalid parameter type!")

}
