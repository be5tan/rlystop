context("Bias variance decomposition")
library(rlystop)

test_that("Bias-variance correctly computes oracle qunatities", {
  # Defining signals
  D <- 10000
  index <- seq(1, D, 1)
  lambda <- index^(-0.5)
  delta <- 0.01
  mu_supersmooth <- 5 * exp(-0.1 * index)
  mu_smooth <- 5000 * abs(sin(0.01 * index)) * index^(-1.6)
  mu_rough <- 250 * abs(sin(0.002 * index))*index^(-0.8)

  # Balanced oracles
  bal_oracle_supersmooth <- bal_oracle(lambda, mu_supersmooth, delta)  # 37
  bal_oracle_smooth      <- bal_oracle(lambda, mu_smooth,      delta)  # 445
  bal_oracle_rough       <- bal_oracle(lambda, mu_rough,       delta)  # 2379

  # 1-balanced oracles
  bal_oracle_supersmooth_1 <- bal_oracle(lambda, mu_supersmooth, delta, alpha = 1)
  # [1] 28
  bal_oracle_smooth_1      <- bal_oracle(lambda, mu_smooth,      delta, alpha = 1)
  # [1] 202
  bal_oracle_rough_1       <- bal_oracle(lambda, mu_rough,       delta, alpha = 1)
  # [1] 706
  
  # Bias-variance decomposition
  B2_supersmooth <- bias2(lambda, mu_supersmooth)
  B2_smooth      <- bias2(lambda, mu_smooth)
  B2_rough       <- bias2(lambda, mu_rough)
  V <- variance(lambda, delta)

  # 1-Bias-variance decomposition
  B2_supersmooth_1 <- bias2(lambda, mu_supersmooth, alpha = 1)
  B2_smooth_1      <- bias2(lambda, mu_smooth,      alpha = 1)
  B2_rough_1       <- bias2(lambda, mu_rough,       alpha = 1)
  V_1 <- variance(lambda, delta, alpha = 1)

  # Check defining properties of the balanced oracle
  m <- bal_oracle_supersmooth
  expect_equal(B2_supersmooth[m - 1] <= V[m - 1], FALSE)
  expect_equal(B2_supersmooth[m] <= V[m], TRUE)
  expect_equal(rbias2(m - 1, lambda, mu_supersmooth) <= rvariance(m - 1, lambda, delta), FALSE)
  expect_equal(rbias2(m, lambda, mu_supersmooth) <= rvariance(m, lambda, delta), TRUE)

  m <- bal_oracle_smooth
  expect_equal(B2_smooth[m - 1] <= V[m - 1], FALSE)
  expect_equal(B2_smooth[m] <= V[m], TRUE)
  expect_equal(rbias2(m - 1, lambda, mu_smooth) <= rvariance(m - 1, lambda, delta), FALSE)
  expect_equal(rbias2(m, lambda, mu_smooth) <= rvariance(m, lambda, delta), TRUE)

  m <- bal_oracle_rough
  expect_equal(B2_rough[m - 1] <= V[m - 1], FALSE)
  expect_equal(B2_rough[m] <= V[m], TRUE)
  expect_equal(rbias2(m - 1, lambda, mu_rough) <= rvariance(m - 1, lambda, delta), FALSE)
  expect_equal(rbias2(m, lambda, mu_rough) <= rvariance(m, lambda, delta), TRUE)

  # Check defining properties of the balanced alpha oracle
  m <- bal_oracle_supersmooth_1
  expect_equal(B2_supersmooth_1[m - 1] <= V_1[m - 1], FALSE)
  expect_equal(B2_supersmooth_1[m] <= V_1[m], TRUE)
  expect_equal(rbias2(m - 1, lambda, mu_supersmooth, 1) <= rvariance(m - 1, lambda, delta, 1), FALSE)
  expect_equal(rbias2(m, lambda, mu_supersmooth, 1) <= rvariance(m, lambda, delta, 1), TRUE)

  m <- bal_oracle_smooth_1
  expect_equal(B2_smooth_1[m - 1] <= V_1[m - 1], FALSE)
  expect_equal(B2_smooth_1[m] <= V_1[m], TRUE)
  expect_equal(rbias2(m - 1, lambda, mu_smooth, 1) <= rvariance(m - 1, lambda, delta, 1), FALSE)
  expect_equal(rbias2(m, lambda, mu_smooth, 1) <= rvariance(m, lambda, delta, 1), TRUE)

  m <- bal_oracle_rough_1
  expect_equal(B2_rough_1[m - 1] <= V_1[m - 1], FALSE)
  expect_equal(B2_rough_1[m] <= V_1[m], TRUE)
  expect_equal(rbias2(m - 1, lambda, mu_rough, 1) <= rvariance(m - 1, lambda, delta, 1), FALSE)
  expect_equal(rbias2(m, lambda, mu_rough, 1) <= rvariance(m, lambda, delta, 1), TRUE)
})
