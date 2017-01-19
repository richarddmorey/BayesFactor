
context('generalTestBF')

data(puzzles)

test_that('generalTestBF works', {
  generalTestBF(RT ~ shape*color*ID, whichRandom="ID", data = puzzles)
})
