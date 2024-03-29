variance <- function(m, lambda, delta, alpha = - 1, filt = c("cutoff", "landw")) {
  #' Variance
  #' 
  #' Computes the variance of the \code{m}-th filter estimator, a given diagonal
  #' design matrix \code{lambda} and noise level delta.
  #'
  #' @param m Integer stopping index.
  #' @param lambda Vector of decreasing, strictly positive entries of the
  #'   diagonal design matrix.
  #' @param delta Numeric noise level.
  #' @param alpha Numeric smoothing parameter. \code{alpha} = -1 gives
  #'   the strong bias. \code{alpha} = 0 gives the weak bias. 
  #' @param filt Character string designating the filter to be used. filt should
  #'   be one of "cutoff" or "landw".
  #'
  #' @return Returns the value of the variancee at index \code{m}.
  #'
  #' @export
  filt <- match.arg(filt)

  if (filt == "cutoff") {
    Vm <- 0 
    if (m > 0) {
      Vm <- sum(lambda[1:m]^(2 * alpha)) * delta^2
    }
  }

  if (filt == "landw") {
    filterTerm <- (1 - (1 - lambda^2)^m)^2
    smoothingTerm <- lambda^(2 * alpha)
    Vm  <- sum(filterTerm * smoothingTerm) * delta^2 
  }

  return(Vm)
}
