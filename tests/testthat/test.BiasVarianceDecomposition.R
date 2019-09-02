context("Bias variance decomposition")
library(rlystop)

test_that("Bias-variance correctly computes oracle quantities", {
  # Defining signals
  D <- 10000
  index <- seq(1, D, 1)
  lambda <- index^(-0.5)
  delta <- 0.01
  muSupersmooth <- 5 * exp(-0.1 * index)
  muSmooth <- 5000 * abs(sin(0.01 * index)) * index^(-1.6)
  muRough <- 250 * abs(sin(0.002 * index))*index^(-0.8)

  # Balanced oracles
  balOracle_supersmooth <- balOracle(lambda, muSupersmooth, delta)  # 37
  balOracle_smooth      <- balOracle(lambda, muSmooth,      delta)  # 445
  balOracle_rough       <- balOracle(lambda, muRough,       delta)  # 2379

  # 1-balanced oracles
  balOracle_supersmooth_1 <- balOracle(lambda, muSupersmooth, delta, 1)  # 28
  balOracle_smooth_1      <- balOracle(lambda, muSmooth,      delta, 1)  # 202
  balOracle_rough_1       <- balOracle(lambda, muRough,       delta, 1)  # 706

  # Bias-variance decomposition
  B2_supersmooth <- sapply(index, bias2,    lambda = lambda, mu = muSupersmooth)
  B2_smooth      <- sapply(index, bias2,    lambda = lambda, mu = muSmooth)
  B2_rough       <- sapply(index, bias2,    lambda = lambda, mu = muRough)
  V              <- sapply(index, variance, lambda = lambda, delta = delta)

  # # 1-Bias-variance decomposition
  B2_supersmooth_1 <- sapply(index, bias2,    lambda = lambda, mu    = muSupersmooth, alpha = 1)
  B2_smooth_1      <- sapply(index, bias2,    lambda = lambda, mu    = muSmooth,      alpha = 1)
  B2_rough_1       <- sapply(index, bias2,    lambda = lambda, mu    = muRough,       alpha = 1)
  V_1              <- sapply(index, variance, lambda = lambda, delta = delta,         alpha = 1)

  # Check defining properties of the balanced oracles
  m <- balOracle_supersmooth
  expect_equal(B2_supersmooth[m - 1] <= V[m - 1], FALSE)
  expect_equal(B2_supersmooth[m] <= V[m], TRUE)

  m <- balOracle_smooth
  expect_equal(B2_smooth[m - 1] <= V[m - 1], FALSE)
  expect_equal(B2_smooth[m] <= V[m], TRUE)

  m <- balOracle_rough
  expect_equal(B2_rough[m - 1] <= V[m - 1], FALSE)
  expect_equal(B2_rough[m] <= V[m], TRUE)

  # Check defining properties of the balanced alpha oracle
  m <- balOracle_supersmooth_1
  expect_equal(B2_supersmooth_1[m - 1] <= V_1[m - 1], FALSE)
  expect_equal(B2_supersmooth_1[m] <= V_1[m], TRUE)

  m <- balOracle_smooth_1
  expect_equal(B2_smooth_1[m - 1] <= V_1[m - 1], FALSE)
  expect_equal(B2_smooth_1[m] <= V_1[m], TRUE)

  m <- balOracle_rough_1
  expect_equal(B2_rough_1[m - 1] <= V_1[m - 1], FALSE)
  expect_equal(B2_rough_1[m] <= V_1[m], TRUE)
})
