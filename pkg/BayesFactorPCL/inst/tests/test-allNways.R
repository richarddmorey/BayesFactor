context("allNways")

test_that("Puzzles example works", {
  data(puzzles)
  
  bf = allNways(y = puzzles$RT, dataFixed = puzzles[,3:4], 
           dataRandom = puzzles$ID, progress = FALSE, extraInfo=TRUE)
  
  expect_that(bf, is_a("data.frame"))
  expect_that(dim(bf), is_equivalent_to(c(5,8)))
})