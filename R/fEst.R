filter_est <- function (m, Y, lambda, m0, filt = c("cutoff")) {
  #' Filter estimator
  #'
  #' Computes the filter estimator with stopping index \code{m} for given data
  #' \code{Y} and diagonal design matrix \code{lambda}.
  #' 
  #' @param m Integer stopping index.
  #' @param Y Numeric vector of observed data.
  #' @param lambda Numeric vector of decreasing, strictly positive entries of
  #'   the diagonal design matrix.
  #' @param m0 Integer starting index.
  #' @param filt Character string giving the filter to be used.
  #'
  #' @return Returns a numeric vector estimating the underlying signal.
  #'
  #' @export
  D <- length(Y)
  mu_hat <- rep(0, D)
  mu_hat[1:m] <- Y[1:m] / lambda[1:m]
  return(mu_hat) 
}
