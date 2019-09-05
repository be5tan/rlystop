context("Classical oracle")
library(rlystop)

test_that("the classical oracle is correct for cutoff", {
  # Defining signals
  D <- 2000
  index <- seq(1, D, 1)
  lambda <- index^(-0.5)
  delta <- 0.01
  muSupersmooth <- 5 * exp(-0.1 * index)
  muSmooth <- 5000 * abs(sin(0.01 * index)) * index^(-1.6)
  muRough <- 250 * abs(sin(0.002 * index))*index^(-0.8)

  filt <- "cutoff"

  # Classical oracles
  claOracle_supersmooth <- claOracle(lambda, muSupersmooth, delta, filt = filt)
  claOracle_smooth      <- claOracle(lambda, muSmooth,      delta, filt = filt)
  claOracle_rough       <- claOracle(lambda, muRough,       delta, filt = filt)

  # 1-classical oracles
  claOracle_supersmooth_1 <- claOracle(lambda, muSupersmooth, delta, 1, filt = filt)  # 28
  claOracle_smooth_1      <- claOracle(lambda, muSmooth,      delta, 1, filt = filt)  # 202
  claOracle_rough_1       <- claOracle(lambda, muRough,       delta, 1, filt = filt)  # 706

  # Bias-variance decomposition
  B2_supersmooth <- sapply(index, bias2,    lambda = lambda, mu = muSupersmooth, filt = filt)
  B2_smooth      <- sapply(index, bias2,    lambda = lambda, mu = muSmooth,      filt = filt)
  B2_rough       <- sapply(index, bias2,    lambda = lambda, mu = muRough,       filt = filt)
  V              <- sapply(index, variance, lambda = lambda, delta = delta,      filt = filt)

  # 1-Bias-variance decomposition
  B2_supersmooth_1 <- sapply(index, bias2,    lambda = lambda, mu    = muSupersmooth, alpha = 1, filt = filt)
  B2_smooth_1      <- sapply(index, bias2,    lambda = lambda, mu    = muSmooth,      alpha = 1, filt = filt)
  B2_rough_1       <- sapply(index, bias2,    lambda = lambda, mu    = muRough,       alpha = 1, filt = filt)
  V_1              <- sapply(index, variance, lambda = lambda, delta = delta,         alpha = 1, filt = filt)

  # Compute the oracles from the decomposition
  testClaOracle_supersmooth <- which.min(B2_supersmooth + V)
  testClaOracle_smooth      <- which.min(B2_smooth + V)
  testClaOracle_rough       <- which.min(B2_rough + V)
  testClaOracle_supersmooth_1 <- which.min(B2_supersmooth_1 + V_1)
  testClaOracle_smooth_1      <- which.min(B2_smooth_1 + V_1)
  testClaOracle_rough_1       <- which.min(B2_rough_1 + V_1)

  # Comparison of the two computations
  expect_equal(testClaOracle_supersmooth,   claOracle_supersmooth)
  expect_equal(testClaOracle_smooth,        claOracle_smooth)
  expect_equal(testClaOracle_rough,         claOracle_rough)
  expect_equal(testClaOracle_supersmooth_1, claOracle_supersmooth_1)
  expect_equal(testClaOracle_smooth_1,      claOracle_smooth_1)
  expect_equal(testClaOracle_rough_1,       claOracle_rough_1)
})

test_that("the classical oracle is correct for Landweber", {
  # Defining signals
  D <- 1000
  index <- seq(1, D, 1)
  lambda <- index^(-0.5)
  delta <- 0.01
  muSupersmooth <- 5 * exp(-0.1 * index)
  muSmooth <- 5000 * abs(sin(0.01 * index)) * index^(-1.6)
  muRough <- 250 * abs(sin(0.002 * index))*index^(-0.8)

  filt <- "landw"

  # Classical oracles
  claOracle_supersmooth <- claOracle(lambda, muSupersmooth, delta, filt = filt)
  claOracle_smooth      <- claOracle(lambda, muSmooth,      delta, filt = filt)
  claOracle_rough       <- claOracle(lambda, muRough,       delta, filt = filt)

  # 1-classical oracles
  claOracle_supersmooth_1 <- claOracle(lambda, muSupersmooth, delta, 1, filt = filt)  # 28
  claOracle_smooth_1      <- claOracle(lambda, muSmooth,      delta, 1, filt = filt)  # 202
  claOracle_rough_1       <- claOracle(lambda, muRough,       delta, 1, filt = filt)  # 706

  # Bias-variance decomposition
  B2_supersmooth <- sapply(index, bias2,    lambda = lambda, mu = muSupersmooth, filt = filt)
  B2_smooth      <- sapply(index, bias2,    lambda = lambda, mu = muSmooth,      filt = filt)
  B2_rough       <- sapply(index, bias2,    lambda = lambda, mu = muRough,       filt = filt)
  V              <- sapply(index, variance, lambda = lambda, delta = delta,      filt = filt)

  # 1-Bias-variance decomposition
  B2_supersmooth_1 <- sapply(index, bias2,    lambda = lambda, mu    = muSupersmooth, alpha = 1, filt = filt)
  B2_smooth_1      <- sapply(index, bias2,    lambda = lambda, mu    = muSmooth,      alpha = 1, filt = filt)
  B2_rough_1       <- sapply(index, bias2,    lambda = lambda, mu    = muRough,       alpha = 1, filt = filt)
  V_1              <- sapply(index, variance, lambda = lambda, delta = delta,         alpha = 1, filt = filt)

  # Compute the oracles from the decomposition
  testClaOracle_supersmooth <- which.min(B2_supersmooth + V)
  testClaOracle_smooth      <- which.min(B2_smooth + V)
  testClaOracle_rough       <- which.min(B2_rough + V)
  testClaOracle_supersmooth_1 <- which.min(B2_supersmooth_1 + V_1)
  testClaOracle_smooth_1      <- which.min(B2_smooth_1 + V_1)
  testClaOracle_rough_1       <- which.min(B2_rough_1 + V_1)

  # Comparison of the two computations
  expect_equal(testClaOracle_supersmooth,   claOracle_supersmooth)
  expect_equal(testClaOracle_smooth,        claOracle_smooth)
  expect_equal(testClaOracle_rough,         claOracle_rough)
  expect_equal(testClaOracle_supersmooth_1, claOracle_supersmooth_1)
  expect_equal(testClaOracle_smooth_1,      claOracle_smooth_1)
  expect_equal(testClaOracle_rough_1,       claOracle_rough_1)
})
