#' Calibration of Computer Models
#'
#' `calibrate` combines simulator and field data and uses given (or default) priors
#' and initial values of the model parameters, runs a full Bayesian Markov Chain Monte Carlo (MCMC)
#' algorithm and samples the posterior distribution of parameters.
#'
#'
#' @details
#'
#' In both matrices of simulation and field data, the first columns are experimental inputs,
#' followed by calibration input columns (only in simulation matrix), and last column is
#' univariate response variable.
#'
#' Posterior distribution of a parameters often do not correspond to well-known
#' distributions. A MCMC algorithm, and especially Metropolis-Hastings (MH)
#' algorithm with adaptive proposal gives a sufficiently close distribution. In
#' addition to point estimates of the parameters, `calibrate` also enables
#' uncertainty quantification for estimated parameters.
#'
#' ## Calibration Model
#'
#' This implementation is based on Kennedy-O'Hagan (KOH) calibration model. In
#' their seminal paper*, Kennedy and O'Hagan augmented simulator output with
#' field observation and fit a three-component model that accounts for simulator
#' input and bias correction using two independent Gaussian Processes (GP) and a
#' third term representing measurement error.
#' \deqn{z_i = \eta (x_i, \kappa) + \delta(x_i) + e_i}
#' In the original paper, the authors proposed a two-stage hierarchical Bayesian
#' model. In the first stage, point estimates of GP hyperparameters are computed
#' using maximum likelihood estimation (MLE) method. In the second stage, these
#' hyperparameters are fixed at their estimated value and run a MCMC algorithm
#' to sample calibration parameters.
#' In contrast, `calibrate` runs MCMC algorithm to sample the posterior
#' distribution of all parameters/hyperparameters. As a result, priors must be
#' specified carefully to reflect the prior expert belief about the distribution
#' and initial values must be chosen as close as possible to prior means.
#'
#' ## Notation
#'
#' parameters of KOH model include all of calibration parameters and will be
#' denoted by \eqn{\kappa}. Since data is scaled and has mean zero, both GPs are
#' specified using their correlation structures and marginal variances. In this
#' implementation, power correlation structure assumed as it has enough
#' flexibility to capture both scale and smoothness of correlation. Therefore,
#' hyperparameters includes scale (\eqn{\theta_S} & \eqn{\theta_B}) and smoothness
#' (\eqn{\alpha_S & \eqn{\alpha_B\eqn{) coefficients and marginal variances (\eqn{\sigma^2_S} &
#' \eqn{\sigma^2_B}) for each GP. Moreover there is a third term that represents
#' measurement error and is specified using its variance (\eqn{\sigma^2_E}).
#'
#'
#' @param sim       A \eqn{m \times (1+p+q)} matrix, representing simulation data,
#' where m, number of rows, is number of simulation runs, p is number of
#' experimental variables, and q is number of calibration variables. Plus one
#' represent the first column, which is the response variable.
#' @param field     A \eqn{n \times (1+p)} matrix, representing field data, where n,
#' number of rows, is number of field observations and p is number of
#' experimental variables.
#' @param Nmcmc     An integer for number of MCMC runs.
#' @param nBurn     An integer for number of MCMC burn ins.
#' @param thining   A double as MCMC thinning rate to remove auto-correlation.
#' @param kappa     A string to specify the prior type of calibration parameters.
#' @param k1        A double as first parameter of the chosen kappa distribution.
#' @param k2        A double as second parameter of the chosen kappa distribution.
#' @param theta     A string to specify the prior type for scale parameters.
#' @param t1        A double as first parameter of the chosen theta distribution.
#' @param t2        A double as second parameter of the chosen theta distribution.
#' @param alpha     A string to specify the prior type for smoothness parameters.
#' @param a1        A double as first parameter of the chosen alpha distribution.
#' @param a2        A double as second parameter of the chosen alpha distribution.
#' @param sigma2    A string to specify the prior type for variance parameters.
#' @param s1        A double as first parameter of the chosen sigma2 distribution.
#' @param s2        A double as second parameter of the chosen sigma2 distribution.
#'
#' @return an output list containing:
#'  * estimates:      Point estimates for all parameters and hyperparameters
#'  * distributions:  \eqn{\frac{(Nmcmc - nBurn)}{thinning} \times k} parameters distribution matrix,
#'                    where k (number of columns) is total number of parameters
#'  * acceptance:     acceptance rate of proposed samples
#'  * logPost:        log posterior distribution of response
#'  * samples:        full unfiltered samples of MCMC runs
#'  * vars:           name of all parameters (based on below notation)
#'  * summary:        summary of all columns of distribution matrix
#'
#' @export
#'
#' @examples examples/calib
#' @references
#' Kennedy MC, O’Hagan A (2001). “Bayesian calibration of computer models.”
#' *Journal of the Royal Statistical Society*, **Series B**, **63(3)**, 425–464
#' <https://www2.stat.duke.edu/~fei/samsi/Oct_09/bayesian_calibration_of_computer_models.pdf>
calibrate <- function(sim, field,
                      Nmcmc  = 11000, nBurn = 1000, thining = 100,
                      kappa  = "gaussian",  k1 = 0.5,  k2 = 0.25,
                      theta  = "logbeta",   t1 = NULL, t2 = NULL,
                      alpha  = "betashift", a1 = 5,    a2 = 1,
                      sigma2 = "jefferys",  s1 = NULL, s2 = NULL) {
  kappaPr    <- setup_prior(kappa,  k1, k2)
  thetaPr    <- setup_prior(theta,  t1, t2)
  alphaPr    <- setup_prior(alpha,  a1, a2)
  sigma2Pr   <- setup_prior(sigma2, s1, s2)
  init       <- setup_cache(sim, field, kappaPr, thetaPr, alphaPr, sigma2Pr)

  mcmc(Nmcmc, nBurn, thining, init, kappaPr, thetaPr, alphaPr, sigma2Pr)

  return(output())
}


