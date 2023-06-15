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
#' parameters of FBC model include all of calibration parameters and will be
#' denoted by \eqn{\kappa}. Since data is scaled and has mean zero, both GPs are
#' specified using their correlation structures and marginal variances. In this
#' implementation, power correlation structure assumed as it has enough
#' flexibility to capture both scale and smoothness of correlation. Therefore,
#' hyperparameters includes scale (\eqn{\theta_S} & \eqn{\theta_B}) and smoothness
#' (\eqn{\alpha_S} & \eqn{\alpha_B}) coefficients and marginal variances (\eqn{\sigma^2_S} &
#' \eqn{\sigma^2_B}) for each GP. Moreover there is a third term that represents
#' measurement error and is specified using its variance (\eqn{\sigma^2_E}).
#'
#'
#' @param sim       A \eqn{m \times (p+q+1)} matrix, representing simulation data,
#' where m, number of rows, is number of simulation runs, p is number of
#' experimental variables, and q is number of calibration variables. Plus one
#' represent the last column, which is the response variable.
#' @param field     A \eqn{n \times (1+p)} matrix, representing field data, where n,
#' number of rows, is number of field observations and p is number of
#' experimental variables.
#' @param Nmcmc     An integer for number of MCMC runs.
#' @param nBurn     An integer for number of MCMC burn ins.
#' @param thining   A double as MCMC thinning rate to remove auto-correlation.
#' @param kappa     A string to specify the prior type of calibration parameters.
#' @param k0        a scaler or vector of same length as number of calibration
#'                  parameters that represent initial value for calibration parameters
#'                  after scaling to (0, 1)
#' @param k1        A double as first parameter of the chosen kappa distribution.
#' @param k2        A double as second parameter of the chosen kappa distribution.
#' @param theta     A string to specify the prior type for scale parameters.
#' @param t0        a scaler in (0, Inf) representing initial values of theta (scale) parameters
#' @param t1        A double as first parameter of the chosen theta distribution.
#' @param t2        A double as second parameter of the chosen theta distribution.
#' @param alpha     A string to specify the prior type for smoothness parameters.
#' @param a0        a scaler in (1, 2) representing initial values of alpha (smoothness) parameters
#' @param a1        A double as first parameter of the chosen alpha distribution.
#' @param a2        A double as second parameter of the chosen alpha distribution.
#' @param sigma2    A string to specify the prior type for variance parameters.
#' @param s0        a scaler or vector of length three that represent initial value for marginal variance
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
calibrate <- function(sim, field,                                           # Data
                      Nmcmc = 2200, nBurn = 200, thining = 20,                  # MCMC
                      kappa   = "beta",         k0  = 0.5, k1 = 1.1,  k2 = 1.1, # Priors
                      thetaS  = "gamma",        ts0 = 0.5, ts1 = 1.1, ts2 = 0.1,
                      alphaS  = "betashift",    as0 = 1.8, as1 = 5,   as2 = 2,
                      thetaB  = "gamma",        tb0 = 0.5, tb1 = 1.1, tb2 = 0.1,
                      alphaB  = "betashift",    ab0 = 1.8, ab1 = 5,   ab2 = 2,
                      sigma2S = "gamma",        ss0 = 1,   ss1 = 0.1, ss2 = 0.1,
                      sigma2B = "inversegamma", sb0 = 1,   sb1 = 0.1, sb2 = 0.1,
                      sigma2E = "inversegamma", se0 = 1,   se1 = 0.1, se2 = 0.1,
                      mu      = "uniform",      m0  = 0,   m1 = -10,  m2 = 10)  {

  init        <- setup_cache(sim, field, Nmcmc,
                             kappa,   k0,  k1,  k2,
                             thetaS,  ts0, ts1, ts2,
                             alphaS,  as0, as1, as2,
                             thetaB,  tb0, tb1, tb2,
                             alphaB,  ab0, ab1, ab2,
                             sigma2S, ss0, ss1, ss2,
                             sigma2B, sb0, sb1, sb2,
                             sigma2E, se0, se1, se2,
                             mu,      m0,  m1,  m2)



  mcmc(Nmcmc, nBurn, thining, init)

  return(output())
}


