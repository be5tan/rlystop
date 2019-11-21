sim <- function (N, lambda, mu, delta, alpha, kappa, filt = 
                 c("cutoff", "landw"), m0 = 0, alphaLoss = - 1) {
  #' Monte Carlo simulation for the filter estimation
  #' 
  #' Performs \code{N} simulation runs of the filter estimation procedure for a
  #' given design matrix \code{lambda}, a given signal \code{mu} and a given
  #' noise level \code{delta} and evaluates the loss and oracle quantities.
  #' 
  #' @param N Integer number of simulation runs.
  #' @param lambda Numeric vector of decreasing, strictly positive entries of
  #'   the diagonal design matrix.
  #' @param mu Numeric vector valued input signal.
  #' @param delta Numeric noise level.
  #' @param alpha Numeric smoothing index for the residuals.
  #' @param kappa Numeric critical stopping value.
  #' @param filt Character string designating the filter to be used. filt should
  #'   be one of "cutoff" or "landw".
  #' @param m0 Integer cut-off index for the two-step procudure.
  #' @param alphaLoss Numeric smoothing parameter for the loss. \code{alpha} =
  #'   -1 gives the strong loss. \code{alpha} = 0 gives the weak loss. 
  #'
  #' @return Returns a data frame containing \code{N} observations of the
  #'  residual stopping rule and loss.
  #' 
  #' @export
  filt <- match.arg(filt)

  D <- length(mu)
  stoppingTime <- rep(0, N)
  loss <- rep(0, N)

  for (iter in 1:N) {
    eps <- stats::rnorm(D, 0, 1)
    Y <- lambda * mu + delta * eps 

    tau <- sTime(Y, lambda, alpha, kappa, filt)
    if (tau < m0) {
      AIC <- rep(0, m0)
      for (m in 1:m0) {
        AIC[m] <- aic(m, Y, lambda, delta, alphaLoss, filt)
      }
      tau2Step <- which.min(AIC)
      stoppingTime[iter] <- tau2Step
      muHat <- fEst(tau2Step, Y, lambda, filt)
    } else {
      muHat <- fEst(tau, Y, lambda, filt)
      stoppingTime[iter] <- tau
    }
    loss[iter] <- sum(lambda^(2 + 2 * alphaLoss) * (muHat - mu)^2)
  }

  df <- data.frame(stoppingTime, loss)
  return(df)
}
