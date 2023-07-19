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
#' field measurement error  variance (\eqn{\sigma^2_E}). If calibration parameters have
#' different initial values, a double vector of same length as number of calibration
#' inputs (determined internally by subtracting the number of columns in `sim` and
#' `field`) can be used.
#'
#' @param sim        \eqn{m \times (1+p+q+1)} numeric matrix, representing simulation
#'                    data, where m is the number of simulation runs, p is number of
#'                    experimental inputs, and q is the number of calibration inputs. The
#'                    plus one column represents the response which is the output of the
#'                    computer code.
#' @param field       \eqn{n \times (1+p)} numeric matrix, representing field data, where n
#'                    is the number of field observations and p is number of experimental
#'                    inputs. Plus one (the first column) represents the field response.
#'
#' @param Nmcmc       integer for number of MCMC runs.
#' @param nBurn       integer for number of MCMC burn ins.
#' @param thinning     integer representing sampling frequency of MCMC results to remove
#'                    auto-correlation.
#'
#'                                 data.
#' @param hypers      nested list containing the priors and initial values for all
#'                    hyperparameters. The notation for the list members are explained
#'                    below. Each member contains four fields (similar to `kappa`): the
#'                    distribution type (`dist`), the initial value (`init`), first
#'                    distribution parameter (`p1`), and second distribution parameter
#'                    (`p2`). The default is a call to `set_hyperPriors()` without
#'                    argument to build the nested list. All arguments of
#'                    `set_hyperPriors()`. In most cases there is no need to change the
#'                    default values. Nevertheless, when needed, user can specify the
#'                    prior distribution for any of the model hyperparameters by supplying
#'                    the changed arguments and their values.
#' @param showProgress  logical indicating whether progress must be displayed at console.
#'                    Default is False.
#' @param kappaDist  string (vector of strings) to specify the prior distribution type(s) for
#'                   calibration parameters.
#' @param kappaInit  double (vector of doubles) that represent initial value(s) of calibration
#'                   parameters to start the MCMC. The default(NA) automates choosing initial value
#'                   for calibration parameters based on their range in simulation
#' @param kappaP1    double (vector of doubles) representing the first parameter(s) of the chosen
#'                   distribution(s)
#' @param kappaP2    double (vector of doubles) representing the first parameter(s) of the chosen
#'                   distribution(s)
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
#' @example man/examples/ex_calibrate.R
#' @references
#' Kennedy MC, O’Hagan A (2001). “Bayesian calibration of computer models.”
#' *Journal of the Royal Statistical Society*, **Series B**, **63(3)**, 425–464
#' <https://www2.stat.duke.edu/~fei/samsi/Oct_09/bayesian_calibration_of_computer_models.pdf>
#' @export
calibrate <- function(sim, field,                                                       # Data
                      Nmcmc  = 2200, nBurn = 200, thinning = 20,                        # MCMC
                      kappaDist = "beta", kappaInit = NA, kappaP1 = 1.1, kappaP2 = 1.1, # Priors
                      hypers = set_hyperPriors(),
                      showProgress = FALSE) {                  # Hyperparameter Priors
  stopifnot((is.matrix(sim) || is.data.frame(sim)),
            (is.matrix(field) || is.data.frame(field)),
            ncol(sim) > 1, ncol(field) > 0)
  stopifnot(Nmcmc > 1, nBurn >= 0, thinning > 0)
  stopifnot(is.list(hypers), names(hypers) == c("thetaS", "alphaS", "thetaB", "alphaB",
                                                "sigma2S", "sigma2B", "sigma2E", "muB"))

  priors <- c(list(kappa = list(dist = kappaDist, init = kappaInit, p1 = kappaP1, p2 = kappaP2)), hypers)
  for (param in names(priors))
    stopifnot(is.character(priors[[param]][['dist']]))
    stopifnot(priors[[param]][['dist']] %in% c("beta", "betashift", "logistic", "gamma",
                                               "uniform", "gaussian", "inversegamma",
                                               "exponential", "jeffreys","lognormal"))
    stopifnot(is.double(priors[[param]][['init']]),
              is.double(priors[[param]][['p1']]),
              is.double(priors[[param]][['p2']]))
    stopifnot(length(priors[[param]][['dist']]) == length(priors[[param]][['init']]))
    stopifnot(length(priors[[param]][['dist']]) == length(priors[[param]][['p1']]))
    stopifnot(length(priors[[param]][['dist']]) == length(priors[[param]][['p2']]))

  init   <- setup_cache(sim, field, priors, Nmcmc, showProgress)
  inds   <- seq(nBurn + 1, Nmcmc, by = thinning)
  output <- mcmc(init, Nmcmc, inds, showProgress)
  return(output)
}


