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
calibrate <- function(sim, field,
                      Nmcmc  = 10000, nBurn = 1, thining = 1,
                      kappa  = "betashift2",  k0 = 0.5,  k1 = 1.2,  k2 = 1.2,
                      theta  = "gamma",       t0 = 0.5,  t1 = 1.5,  t2 = 0.1,
                      alpha  = "betashift",   a0 = 1.75, a1 = 5,    a2 = 2,
                      sigma2 = "gamma",       s0 = 1,    s1 = 1.5,  s2 = 1.5) {
  args    <- formals(calibrate)
  Nmcmc   <- if (!is.null(args$Nmcmc))   Nmcmc   else 20000
  nBurn   <- if (!is.null(args$nBurn))   nBurn   else 1
  thining <- if (!is.null(args$thining)) thining else 1
  kappa   <- if (!is.null(args$kappa))   kappa   else "betashift2"
  k1      <- if (!is.null(args$k1))      k1      else 1.2
  k2      <- if (!is.null(args$k2))      k2      else 1.2
  theta   <- if (!is.null(args$theta))   theta   else "gamma"
  t1      <- if (!is.null(args$t1))      t1      else 1.5
  t2      <- if (!is.null(args$t2))      t2      else 0.1
  alpha   <- if (!is.null(args$alpha))   alpha   else "betashift"
  a1      <- if (!is.null(args$a1))      a1      else 5
  a2      <- if (!is.null(args$a2))      a2      else 2
  sigma2  <- if (!is.null(args$sigma2))  sigma2  else "gamma"
  s1      <- if (!is.null(args$s1))      s1      else 1.5
  s2      <- if (!is.null(args$s2))      s2      else 1.5

  kappaPr          <- setup_prior(kappa,  k1, k2)
  thetaPr          <- setup_prior(theta,  t1, t2)
  alphaPr          <- setup_prior(alpha,  a1, a2)
  sigma2Pr         <- setup_prior(sigma2, s1, s2)
  init             <- setup_cache(sim, field, kappaPr, thetaPr, alphaPr, sigma2Pr, k0, t0, a0, s0)
  cache$Phi        <- matrix(nrow = Nmcmc, ncol = cache$k)
  cache$Phi[1, ]   <- init$phi
  cache$logPost    <- double(Nmcmc)
  cache$logPost[1] <- init$logPost

  cache$priors <- list(kappa  = c(kappa,  k1, k2),
                       theta  = c(theta,  t1, t2),
                       alpha  = c(alpha,  a1, a2),
                       sigma2 = c(sigma2, s1, s2))

  cat("initial values: ", round(init$phi, 3), "\n")

  mcmc(Nmcmc, nBurn, thining, kappaPr, thetaPr, alphaPr, sigma2Pr)

  return(output())
}


