aic <- function(m, Y, lambda, delta, alphaLoss = - 1, filt = c("cutoff", "landw")) {
  #' Akaike information criterion
  #' 
  #' Computes the AIC of the \code{m}-th filter estimator.
  #'
  #' @param m Integer stopping index.
  #' @param Y Numeric vector of observed data.
  #' @param lambda Vector of decreasing, strictly positive entries of the
  #'   diagonal design matrix.
  #' @param delta Numeric noise level.
  #' @param alphaLoss Numeric smoothing parameter for the loss. \code{alpha} =
  #'   -1 gives the strong loss. \code{alpha} = 0 gives the weak loss. 
  #' @param filt Character string designating the filter to be used. filt should
  #'   be one of "cutoff" or "landw".
  #'
  #' @return Returns the value of the AIC at index \code{m}.
  #'
  #' @export
  filt <- match.arg(filt)

  muHat <- fEst(m, Y, lambda, filt)
  # rss <- sum(lambda^(2 + 2 * alphaLoss) * (Y - lambda * muHat)^2)
  rss <- - sum(lambda[1:m]^(2 * alphaLoss) * Y[1:m]^2)
  aic <- rss + 2 * variance(m, lambda, delta, alphaLoss, filt)
  return(aic)
}
