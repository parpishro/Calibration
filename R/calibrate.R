#' Calibration of Computer Models
#'
#' `calibrate` uses a Bayesian framework to estimate the posterior density distribution of
#' calibration model parameters. It takes in simulator data, and field data, and parameter
#' prior specifications as arguments, and forms a numeric joint distribution of all parameters as
#' output.
#'
#'
#' @details
#'
#' Both simulation and field data are required to be dataframe or matrix. The first column in both
#' datasets is the response. In field data rest of columns are experimental inputs and in simulation
#' data there are additional columns at the end that represent calibration inputs.
#'
#' Posterior distribution of a parameters often do not correspond to a well-known distributions. A
#' version of MCMC algorithm, called Metropolis-Within-Gibbs algorithm, with adaptive proposal is
#' used to draw from joint posterior distribution of parameters. The full Bayesian approach of
#' `calibrate` provides a measure of uncertainty quantification as well as point estimates for the
#' parameters. More information about calibration model and package functionality can be found in
#' package vignette.
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
#' @param nMCMC       integer for number of MCMC runs.
#' @param nBurn       integer for number of MCMC burn ins.
#' @param thinning     integer representing sampling frequency of MCMC results to remove
#'                    auto-correlation.
#'
#'                                 data.
#' @param hypers      nested list containing the priors and initial values for all
#'                    hyperparameters. The notation for the list members are explained
#'                    below. Each member contains four fields (similar to `kappa`): the
#'                    distribution type (`dist`), first
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
#'  * priorFamilies:  prior specifications of all parameters
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
calibrate <- function(sim, field,                                      # Data
                      nMCMC  = 1100, nBurn = 100, thinning = 10,       # MCMC
                      kappaDist = "beta", kappaP1 = 1.1, kappaP2 = 1.1, # Calibration Prior(s)
                      hypers = set_hyperPriors(),                      # Hyperparameter Priors
                      showProgress = FALSE) {
  if (!is.matrix(sim) && !is.data.frame(sim))
     stop("'sim' argument must be either matrix or data frame")
  if (!is.matrix(field) && !is.data.frame(field))
     stop("'field' argument must be either matrix or data frame")
  if (ncol(sim) < 2 || ncol(field) < 1)
     stop("'sim' and `field` arguments must have at least 2 & 1 columns respectively!")
  if (nMCMC <= 1)
     stop("'nMCMC' must be a integer and larger than 1!")
  if (nBurn < 0)
     stop("'nBurn' must be a non-negative integer!")
  if (thinning <= 0)
     stop("'thinning' must be a positive integer!")
  if (!is.list(hypers) || sum(names(hypers) != c("thetaS", "alphaS", "thetaB", "alphaB", "sigma2S",
                                             "sigma2B", "sigma2E", "muB")) != 0)
     stop("'hyper' argument does not have the proper format. Use `set_hyperPriors()` function to
          set hyperparameter priors!")


  priorFamilies <- c(list(kappa = list(dist = kappaDist, p1 = kappaP1, p2 = kappaP2)), hypers)
  for (fam in priorFamilies) {
    if (!is.character(fam$dist))
      stop("'dist' argument must be character string!")
    if (sum(!(fam$dist %in% c("uniform", "normal", "normalTr", "lognormal", "gamma", "inversegamma",
                              "beta", "betashift", "logbeta", "logistic","exponential", "fixed"))) > 0)
      stop("'dist' argument must be one of 'uniform', 'normal', 'normalTr', 'lognormal', 'gamma',
           'inversegamma', 'beta', 'betashift', 'logbeta', 'logistic','exponential', 'fixed'")
    if (sum(!is.double(fam$p1)) + sum(!is.double(fam$p2)) > 0)
      stop("Both 'p1' and `p2` arguments must be double!")
    if (length(fam$dist) != length(fam$p1) || length(fam$dist) != length(fam$p2))
      stop("Length of 'dist', 'p1', and 'p2' arguments must be same!")
  }

  init   <- setup_cache(sim, field, priorFamilies, nMCMC, showProgress)
  inds   <- seq(nBurn + 1, nMCMC, by = thinning)
  output <- mcmc(init, nMCMC, inds, showProgress)
  return(output)
}


