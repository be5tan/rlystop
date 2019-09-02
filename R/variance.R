variance <- function(m, lambda, delta, alpha = - 1, filt = c("cutoff")) {
  #' Variance
  #' 
  #' Computes the variance of the filter estimator a given diagonal design
  #' matrix \code{lambda} and noise level delta.
  #'
  #' @param m Integer stopping index.
  #' @param lambda Vector of decreasing, strictly positive entries of the
  #'   diagonal design matrix.
  #' @param delta Numeric noise level.
  #' @param alpha Numeric smoothing parameter.
  #' @param filt Character string giving the filter to be used.
  #'
  #' @return Returns all values of the variancee in a double vector.
  #'
  #' @export
  if (filt == "cutoff") {
    Vm <- sum(lambda[1:m]^(2 * alpha)) * delta^2
  }
  if (filt == "landw") {
    filterTerm <- (1 - (1 - lambda^2)^m)^2
    smoothingTerm <- lambda^(2 * alpha)
    Vm  <- sum(filterTerm * smoothingTerm) * delta^2 
  }
}
