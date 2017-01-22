
context('generalTestBF')

data(puzzles)

test_that('generalTestBF works', {
  bf <- generalTestBF(RT ~ shape*color*ID, whichRandom="ID", data = puzzles)
  expect_that(bf, is_a("BFBayesFactor"))
  expect_that(length(bf), is_equivalent_to(18))
})
