
context("contingencyBF")

test_that("contingencyBF works", {
  data <- matrix(c(3, 6, 4, 9), nrow=2)
  contingencyTableBF(data, sampleType='poisson')
  contingencyTableBF(data, sampleType='jointMulti', fixedMargin='rows')
  contingencyTableBF(data, sampleType='jointMulti', fixedMargin='cols')
  contingencyTableBF(data, sampleType='indepMulti', fixedMargin='rows')
  contingencyTableBF(data, sampleType='indepMulti', fixedMargin='cols')
  contingencyTableBF(data, sampleType='hypergeom')
})
