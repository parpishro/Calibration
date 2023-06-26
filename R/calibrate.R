#' Calibration of Computer Models
#'
#' `calibrate` uses a Bayesian framework to estimate the posterior density distribution of
#' calibration model parameters. It takes in simulator data, and field data, and parameter
#' priors as arguments, and forms a numeric joint distribution of all parameters as output.
#'
#'
#' @details
#'
#' Both simulation and field data are required to be dataframe or matrix. The first column
#' in both datasets is the response. In field data rest of columns are experimental inputs
#' and in simulation data there are additional columns at the end that represent
#' calibration inputs.
#'
#' Posterior distribution of a parameters often do not correspond to well-known
#' distributions. A MCMC algorithm, and especially Metropolis-Hastings (MH)
#' algorithm with adaptive proposal is used to draw from joint posterior distribution of
#' parameters. The full Bayesian approach of `calibrate` enables both point estimates and
#' inference about the uncertainty in parameter estimates. More information about
#' calibration model and package functionality can be found in package vignette.
#'
#'
#' ## Notation
#'
#' parameters of FBC model include all of calibration parameters and will be
#' denoted by \eqn{\kappa}. Both simulator and bias-correction processes are modeled
#' using Gaussian Processes (GP), which in turn are specified with their respective
#' correlation structures and marginal variances. In this implementation, power
#' correlation structure is assumed for its flexibility. In summary, model parameters are
#' calibration parameters (\eqn{\kappa}), simulator GP hyperparameters (marginal variance
#' \eqn{\sigma^2_S}, scale \eqn{\theta_S}, and smoothness \eqn{\alpha_S}), bias-correction
#' GP hyperparameters (\eqn{\sigma^2_B}, \eqn{\theta_B}, and \eqn{\alpha_B}), and finally
#' field measurement error  variance (\eqn{\sigma^2_E}). If
#'                  calibration parameters have different initial values, a double vector
#'                  of same length as number of calibration inputs (determined internally
#'                  by subtracting the number of columns in `sim` and `field`) can be used.
#'
#' @section Data:
#'
#' @param sim     \eqn{m \times (1+p+q+1)} numeric matrix, representing simulation data,
#' where m is the number of simulation runs, p is number of experimental inputs, and q is
#' thenumber of calibration inputs. The plus one column represents the response which is
#' the output of the computer code.
#' @param field  \eqn{n \times (1+p)} numeric matrix, representing field data, where n is
#' the number of field observations and p is number of experimental inputs. Plus one (the
#' first column) represents the field response.
#'
#' @section MCMC:
#'
#' @param Nmcmc     integer for number of MCMC runs.
#' @param nBurn     integer for number of MCMC burn ins.
#' @param thining   integer representing sampling frequency of MCMC results to remove
#'                  auto-correlation.
#'
#' @section Priors:
#'
#' @param kappa     string (vector of strings) to specify the prior distribution type(s)
#'                  for calibration parameters.
#' @param k0        double (vector of doubles) representing initial value(s) for
#' @param p1        double (vector of doubles) representing the first parameter(s) of the
#'                  chosen distribution(s)
#' @param p2        double (vector of doubles) representing the first parameter(s) of the
#'                  chosen distribution(s)
#' @param hypers    a nested list containing the priors and initial values for all
#'                  hyperparameters. The notation for the list members are explained below.
#'                  Each member contains four fields: the distribution type (`dist`), the
#'                  initial value (`init`), first distribution parameter (`p1`), and
#'                  second distribution parameter (`p2`). The default list (`hyperPriors`)
#'                  fully specifies the priors for all hyperparameters of the calibration
#'                  model except calibration parameter(s). They are specified based on
#'                  literature consensus and in most cases user does not need to change
#'                  them. However, when needed, user can change the priors using
#'                  `control()` function. Note that the default value will be changed until next
#'                  restart of the package.
#'
#'
#' @return a `fbc` object containing:
#'  * Phi:            A numeric matrix in which each row represents a draw from joint
#'                    posterior distribution of parameters (after thinning to remove
#'                    autocorrelation between consecutive draws) and each column
#'                    represents a parameter of the model. In essence each column
#'                    approximates the marginal distribution of that parameter.
#'  * estimates:      a dataframe containing the summary of estimates of parameters and
#'                    their uncertainty
#'  * logPost:        a numeric vector of same length as Phi rows representing the log of
#'                    posterior likelihood
#'  * priors:         prior specifications of all parameters
#'  * acceptance:     acceptance rate of MH MCMC algorithm
#'  * vars:           name of all parameters (based on below notation)
#'  * cache:          an environment containing original datasets and indexing variables
#'                    that is used in `predict()` function.
#'
#' @export
#'
#' @examples examples/calib (TODO)
#' @references
#' Kennedy MC, O’Hagan A (2001). “Bayesian calibration of computer models.”
#' *Journal of the Royal Statistical Society*, **Series B**, **63(3)**, 425–464
#' <https://www2.stat.duke.edu/~fei/samsi/Oct_09/bayesian_calibration_of_computer_models.pdf>
calibrate <- function(sim, field,                                                      # Data
                      Nmcmc=1100, nBurn=100, thinning=1,                              # MCMC
                      kappa="beta", k0=NA, p1=1.1, p2=1.1, hypers=set_hyperPriors()) { # Priors
  ical  <- (ncol(field)+1):ncol(sim)
  if (is.na(k0))
    k0   <- apply(sim[, ical, drop=F], 2, mean)
  priors <- c(list(kappa=list(dist=kappa, init=k0, p1=p1, p2=p2)), hypers)
  init   <- setup_cache(sim, field, priors, Nmcmc)
  inds   <- seq(nBurn+1, Nmcmc, by=thinning)
  output <- mcmc(init, Nmcmc, inds)
  return(output)
}


