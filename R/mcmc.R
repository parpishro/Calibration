#' Adaptive Metropolis-within-Gibbs MCMC Sampling
#'
#' Runs an adaptive Metropolis-within-Gibbs MCMC that computes posterior log likelihood in
#' each run and updates the model parameters. New parameters are proposed according to a
#' symmetric and adaptive proposal function. After removing burn-in runs and thinning the
#' samples, the columns of `Phi` matrix in the output approximates the marginal posterior
#' distribution of parameters.
#'
#' @param init     list containing initialized `Phi` mtrix and `logPost` vector
#' @param nMCMC    integer representing total number of MCMC runs
#' @param inds     indices of final MCMC draws after burn-in and thinning
#' @param showProgress  logical indicating whether progress must be displayed at console.
#'                    Default is False.
#'
#' @importFrom stats runif
#' @noRd
#' @return fbc object containing:
#'            * Phi         matrix of posterior estimates: each row represents a draw from
#'                          joint posterior distribution and each column represent
#'                          estimates for marginal distribution of a model parameter
#'            * estimates   summary of estimates for all parameters and their uncertainty
#'            * logPost     vector of log posterior likelihood values
#'            * priors      a nested list representing prior specification
#'            * acceptance  vector representing acceptance rate of each parameter
#'            * vars        name of parameters (column headers of `Phi`)
#'            * cache       environment containing original data and indices
mcmc <- function(init, nMCMC, inds, showProgress) {
  Phi        <- init$Phi
  logPost    <- init$logPost
  priorFns   <- cache$priorFns
  l          <- cache$l
  ifixed     <- cache$ifixed
  accepLast  <- double(l)
  accepted   <- double(l)
  accepRate  <- double(l)
  sdProp     <- double(l)
  res        <- cache$y - Phi[1, cache$imuB]
  for (i in 2:nMCMC) {
    logPost[i]  <- logPost[i - 1]
    for (j in 1:l) {
      if (j %in% ifixed) next
      out       <- proposal(Phi[i - 1, j], j, i, accepRate[j], sdProp[j])
      changed   <- out$proposed
      sdProp[j] <- out$sd

      if (j == 1) {
        params <- c(changed, Phi[i - 1, 2:l])
        update_cov(params, ichanged = j)
      } else if (j == l) {
        params <- c(Phi[i, 1:(l - 1)], changed)
        res    <- cache$y - params[cache$imuB]
      } else {
        params <- c(Phi[i, 1:(j - 1)], changed, Phi[(i - 1), (j + 1):l])
        update_cov(params, ichanged = j)
      }

      logPrior       <- 0
      for (h in 1:l)
        logPrior     <- logPrior + priorFns[[h]](params[h])
      lPost  <- logPrior - 0.5*(cache$logDetCov + drop(t(res) %*% (cache$InvCov) %*% res))


      if (is.finite(lPost) && (lPost - logPost[i] > log(runif(1)))) {
        Phi[i, j]   <- changed
        logPost[i]  <- lPost
        accepted[j] <- accepted[j] + 1
      } else {
        Phi[i, j]      <- Phi[i - 1, j]
      }
    }

    if (i == 5) {
      accepRate  <- (accepted - accepLast)/4
      accepLast  <- accepted
    }
    if (i > 5 && i %% 10 == 5) {
      accepRate  <- (accepted - accepLast)/10
      accepLast  <- accepted
    }


    if (showProgress && i %% floor(nMCMC/100) == 0 && i/floor(nMCMC/100) <= 100) {   #
      cat("------------------------------------------------------------------------", "\n")
      cat("finished ",  (i/floor(nMCMC/100)), "% of MCMC runs...", "\n")
      cat("acceptance ratio in the last batch:   ", round(accepRate, 2), "\n")
      cat("parameter sample:                     ", round(Phi[i, ], 2), "\n")
      cat("proposal sd:                          ", round(sdProp, 2), "\n")
      cat("total acceptance ratio:               ", round(accepted/i, 2), "\n")
    }
  }
  if (showProgress) {
    cat("------------------------------------------------------------------------", "\n")
    cat("Completed MCMC runs.", "\n")
    cat("total acceptance ratio ", round(accepted/nMCMC, 2), "\n")
  }

  cache$Params     <- Phi[inds, ]
  cache$logPost    <- logPost[inds]
  cache$acceptance <- accepted/nMCMC
  return(output())

}
