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
  #' @param alphaLoss Numeric smoothing parameter for the loss.
  #' @param filt Character string giving the filter to be used.
  #'
  #' @return Returns the value of the squared bias at index \code{m}.
  #'
  #' @export
  filt <- match.arg(filt)

  muHat <- fEst(m, Y, lambda, filt)
  rss <- sum(lambda^(2 + 2 * alphaLoss) * (Y - lambda * muHat)^2)
  aic <- rss + 2 * variance(m, lambda, delta, alphaLoss, filt)
}
