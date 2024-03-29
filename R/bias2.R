bias2 <- function(m, lambda, mu, alpha = - 1, filt = c("cutoff", "landw")) {
  #' Squared bias
  #' 
  #' Computes the squared bias of the \code{m}-th filter estimator, a given
  #' diagonal design matrix \code{lambda} and signal \code{mu}.
  #'
  #' @param m Integer stopping index.
  #' @param lambda Vector of decreasing, strictly positive entries of the
  #'   diagonal design matrix.
  #' @param mu Vector valued input signal. 
  #' @param alpha Numeric smoothing parameter. \code{alpha} = -1 gives
  #'   the strong bias. \code{alpha} = 0 gives the weak bias. 
  #' @param filt Character string designating the filter to be used. filt should
  #'   be one of "cutoff" or "landw".
  #'
  #' @return Returns the value of the squared bias at index \code{m}.
  #'
  #' @export
  filt <- match.arg(filt)

  if (filt == "cutoff") {
    D <- length(mu) 
    B2m <- sum(lambda[(m + 1):D]^(2 + 2 * alpha) * mu[(m + 1):D]^2)
  }

  if (filt == "landw") {
    filterTerm <- (1 - lambda^2)^(2 * m)
    smoothingTerm <- lambda^(2 + 2 * alpha)
    B2m <- sum(filterTerm * smoothingTerm * mu^2)
  }

  return(B2m)
}
