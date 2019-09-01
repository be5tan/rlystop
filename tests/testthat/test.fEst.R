context("Filter estimators")
library(rlystop)

test_that("elementary cuf-off estimation works", {
  D <- 10
  lambda <- seq(1, D)
  mu <- rnorm(D, 0, 1)
  Y <- lambda * mu
  for (m in 1:D) {
    expect_equal(fEst(m, Y, lambda), c(mu[1:m], rep(0, D - m)))
  }
})

