sim <- function (N, lambda, mu, delta, alpha, kappa, filt = c("cutoff"), m0 = 0,
                 alphaLoss = - 1) {
  #' Monte Carlo simulation for the filter estimation
  #' 
  #' Performs \code{num_run} simulation runs of the filter estimation
  #' procedure for a given design matrix \code{lambda}, a given signal
  #' \code{mu} and a given noise level \code{delta} and evaluates the loss and
  #' oracle quantities.
  #' 
  #' @param N Integer number of simulation runs.
  #' @param lambda Numeric vector of decreasing, strictly positive entries of
  #'   the diagonal design matrix.
  #' @param mu Numeric vector valued input signal.
  #' @param delta Numeric positive level of random noise.
  #' @param alpha Numeric smoothing index for the residuals.
  #' @param kappa Numeric positive critical stopping value.
  #' @param filt Character string giving the filter to be used.
  #' @param m0 Integer cut-off index for the two-step procudure.
  #' @param alphaLoss Numeric smoothing index for the loss. 
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
    stoppingTime[iter] <- tau
    muHat <- fEst(tau, Y, lambda, filt)
    loss[iter] <- sum(lambda^(1 + alphaLoss) * (muHat - mu)^2)
  }
  df <- data.frame(stoppingTime, loss)

  return(df)
}
