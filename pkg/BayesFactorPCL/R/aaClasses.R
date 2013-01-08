# https://stat.ethz.ch/pipermail/r-devel/2010-May/057506.html
## for 'i' in x[i] or A[i,] : (numeric = {double, integer})
setClassUnion("index", members =  c("numeric", "logical", "character"))

#setClassUnion("missOrLog", members =  c("missing", "logical"))


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

setOldClass("mcmc")
setClass("BFmcmc", contains = "mcmc", representation(model="BFmodel",data = "data.frame"))
