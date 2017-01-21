
context("correlationBF")

set.seed(0)

test_that("correlation works", {
  x <- runif(100)
  y <- x + rnorm(100)
  correlationBF(x, y)
})
