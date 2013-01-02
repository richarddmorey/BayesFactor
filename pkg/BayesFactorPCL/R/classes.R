
setClass("BFmodel", representation(
  type = "character",
  identifier = "list",
  prior = "list",
  dataTypes = "character",
  shortName = "character",
  longName = "character",
  version = "character"
))

setClass("BFBayesFactor", representation(
  numerator = "list",
  denominator = "BFmodel",
  bayesFactor = "data.frame",
  data = "data.frame",
  version = "character"
))

setClass("BFlinearModel", contains = "BFmodel")
setClass("BFoneSample", contains = "BFlinearModel")
setClass("BFindepSample", contains = "BFlinearModel")

setClass("BFBayesFactorList", contains = "list", representation(version="character"))