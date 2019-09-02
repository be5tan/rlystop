bias2 <- function(m, lambda, mu, alpha = - 1, filt = c("cutoff")) {
  #' Squared bias
  #' 
  #' Computes the squared bias of the filter estimator a given diagonal design
  #' matrix \code{lambda} and signal ' \code{mu}.
  #'
  #' @param m Integer stopping index.
  #' @param lambda Vector of decreasing, strictly positive entries of the
  #'   diagonal design matrix.
  #' @param mu Vector valued input signal. 
  #' @param alpha Numeric smoothing parameter.
  #' @param filt Character string giving the filter to be used.
  #'
  #' @return Returns the value of the squared bias at index \code{m}.
  #'
  #' @export
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
