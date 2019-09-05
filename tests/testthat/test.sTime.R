library(rlystop)
context("Residual stopping time")

test_that("Smoothed residual stopping time is correct for cut-off", {
  filt = "cutoff"
  D <- 100
  index <- seq(1, D, 1)
  lambda <- index^(-0.5)
  mu <- 250 * abs(sin(0.002 * index)) * index^(-0.8)
  Y <- lambda * mu 
  kappa <- 1
  for (alpha in c(0, 0.5, 1)) {
    for (m in 1:D) {
      mu_hat <- fEst(m, Y, lambda, filt)  
      residuals2 <- sum(lambda^(2 * alpha) * (Y - lambda * mu_hat)^2)
      expect_equal(sTime(Y, lambda, alpha, kappa, filt) <= m, residuals2 <=
                   kappa)
    }
  }

  # High dimensional edge case
  delta <- 1e-05
  D <- 1e06
  index <- seq(1, D, 1)
  lambda <- index^(-0.5)
  mu <- 5 * exp(-0.1 * index)
  eps <- rnorm(D, 0, 1)
  Y <- lambda * mu + delta * eps
  alpha <- 0
  kappa <- sum(lambda^(2 * alpha)) * delta^2 
  m <- 1
  residuals2 <- sum(lambda[m:D]^(2 * alpha) * Y[m:D]^2)
  while(residuals2 > kappa) {
    residuals2  <- residuals2 - lambda[m]^(2 * alpha) * Y[m]^2 
    m <- m + 1 
  }
  expect_equal(sTime(Y, lambda, alpha, kappa, filt), m - 1)
})

test_that("Smoothed residual stopping time is correct for Landweber", {
  filt = "landw"
  D <- 100
  index <- seq(1, D, 1)
  lambda <- index^(-0.5)
  mu <- 250 * abs(sin(0.002 * index)) * index^(-0.8)
  Y <- lambda * mu 
  kappa <- 1
  for (alpha in c(0, 0.5, 1)) {
    for (m in 1:D) {
      mu_hat <- fEst(m, Y, lambda, filt)  
      residuals2 <- sum(lambda^(2 * alpha) * (Y - lambda * mu_hat)^2)
      expect_equal(sTime(Y, lambda, alpha, kappa, filt) <= m, residuals2 <=
                   kappa)
    }
  }
})
