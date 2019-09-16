fEst <- function (m, Y, lambda, filt = c("cutoff", "landw")) {
  #' Filter estimator
  #'
  #' Computes the filter estimator with stopping index \code{m} for given data
  #' \code{Y} and diagonal design matrix \code{lambda}.
  #' 
  #' @param m Integer stopping index.
  #' @param Y Numeric vector of observed data.
  #' @param lambda Numeric vector of decreasing, strictly positive entries of
  #'   the diagonal design matrix.
  #' @param filt Character string giving the filter to be used.
  #'
  #' @return Returns a numeric vector estimating the underlying signal.
  #'
  #' @export
  filt <- match.arg(filt)

  if (filt == "cutoff") {
    D <- length(Y)
    muHat <- rep(0, D)
    muHat[1:m] <- Y[1:m] / lambda[1:m]
  }

  if(filt == "landw") {
    muHat <- landw(m, Y, lambda)
  } 

  return(muHat)
}
