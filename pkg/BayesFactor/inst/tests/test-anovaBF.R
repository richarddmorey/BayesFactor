
context("anovaBF")

test_that("Puzzles example works", {
  data(puzzles)

  bf = anovaBF(RT ~ shape*color + ID, data = puzzles, whichRandom = "ID",
               iterations=1000, progress = FALSE)

  expect_that(bf, is_a("BFBayesFactor"))
  expect_that(length(bf), is_equivalent_to(4))
})
