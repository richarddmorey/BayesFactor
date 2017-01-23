
context("regressionBF")

set.seed(0)

test_that("regressionBF works", {
  x <- rnorm(100)
  a <- rnorm(100)
  y <- x + a + rnorm(100)
  data <- data.frame(y, x, a)
  regressionBF(y ~ x + a, data)
})

test_that("regressionBF errors appropriately", {
  x <- rnorm(100)
  a <- rnorm(100)
  y <- x + a + rnorm(100)
  data <- data.frame(y, x, a)
  expect_error(
    regressionBF(y ~ x * a, data),
    'Interactions not allowed in regressionBF (try generalTestBF)',
    fixed=TRUE)
})
