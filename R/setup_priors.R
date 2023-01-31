setup_priors <- function(t, o, a, s) {
  if (is.na(t) || is.na(o) || is.na(a) || is.na(s))
    stop("at least one prior distribution is missing!")
  if (t == -1 || o == -1 || a == -1 || s == -1)
    stop("at least one prior distribution is invalid")
  if (t > 2 || o > 2 || a > 2 || s > 2)
    stop("at least one prior distribution is still under development")
  return(TRUE)
}
