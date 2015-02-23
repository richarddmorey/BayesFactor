# https://stat.ethz.ch/pipermail/r-devel/2010-May/057506.html
## for 'i' in x[i] or A[i,] : (numeric = {double, integer})
# setClassUnion("index", members =  c("numeric", "logical", "character"))



#' General S4 classes for representing models for comparison
#' 
#' The \code{BFmodel} is a general S4 class for representing models for comparison. The more classes 
#' \code{BFlinearModel}, \code{BFindepSample}, and \code{BFoneSample} inherit directly from \code{BFmodel}. 
#' 
#'   \describe{
#'   These model classes all have the following slots defined:
#'    \item{type}{Model type}
#'    \item{identifier}{a list uniquely identifying the model from other models of the same type}
#'    \item{prior}{list giving appropriate prior settings for the model}
#'    \item{dataTypes}{a character vector whose names are possible columns in the data; elements specify the corresponding data type, currently one of c("fixed","random","continuous")}
#'    \item{shortName}{a short, readable identifying string}
#'    \item{longName}{a longer, readable identifying string}
#'    \item{analysis}{object storing information about a previous analysis of this model}
#'    \item{version}{character string giving the version and revision number of the package that the model was created in}
#'    }
#' @name BFmodel-class
#' @rdname model-classes
#' @export
setClass("BFmodel", representation(
  type = "character",
  identifier = "list",
  prior = "list",
  dataTypes = "character",
  shortName = "character",
  longName = "character",
  analysis = "list",
  version = "character"
))

#' @name BFcontingencyTable-class
#' @rdname model-classes
setClass("BFproportion", contains = "BFmodel")

#' @name BFcontingencyTable-class
#' @rdname model-classes
setClass("BFcontingencyTable", contains = "BFmodel")

#' @name BFlinearModel-class
#' @rdname model-classes
setClass("BFlinearModel", contains = "BFmodel")

#' @name BFoneSample-class
#' @rdname model-classes
setClass("BFoneSample", contains = "BFlinearModel")

#' @name BFoneSample-class
#' @rdname model-classes
setClass("BFmetat", contains = "BFmodel")

#' @name BFindepSample-class
#' @rdname model-classes
setClass("BFindepSample", contains = "BFlinearModel")

#' General S4 class for representing multiple Bayes factor model comparisons, all against the same model
#' 
#' The \code{BFBayesFactor} class is a general S4 class for representing models model comparison via Bayes factor.
#'   
#'   \code{BFBayesFactor} objects can be inverted by taking the reciprocal and can
#'    be divided by one another, provided both objects have the same denominator. In addition, 
#'    the \code{t} (transpose) method can be used to invert Bayes factor objects.

#'   \describe{
#'   The \code{BFBayesFactor} class has the following slots defined:
#'    \item{numerator}{a list of models all inheriting \code{BFmodel}, each providing a single denominator}
#'    \item{denominator}{a single \code{BFmodel} object serving as the denominator for all model comparisons}
#'    \item{bayesFactor}{a data frame containing information about the comparison between each numerator and the denominator}
#'    \item{data}{a data frame containing the data used for the comparison}
#'    \item{version}{character string giving the version and revision number of the package that the model was created in}
#'    }
#' @name BFBayesFactor-class
#' @export
#' @examples
#' ## Compute some Bayes factors to demonstrate division and indexing
#' data(puzzles)
#' bfs <- anovaBF(RT ~ shape*color + ID, data = puzzles, whichRandom = "ID", progress=FALSE)
#' 
#' ## First and second models can be separated; they remain BFBayesFactor objects
#' b1 = bfs[1]
#' b2 = bfs[2]
#' b1
#' 
#' ## We can invert them, or divide them to obtain new model comparisons
#' 1/b1
#' b1 / b2
#' 
#' ## Use transpose to create a BFBayesFactorList
#' t(bfs) 
setClass("BFBayesFactor", representation(
  numerator = "list",
  denominator = "BFmodel",
  bayesFactor = "data.frame",
  data = "data.frame",
  version = "character"
))

#' General S4 class for representing a collection of Bayes factor model 
#' comprisons, each against a different denominator
#' 
#' The \code{BFBayesFactorList} class is a general S4 class for representing 
#' models model comparison via Bayes factor. See the examples for demonstrations
#' of BFBayesFactorList methods.
#' 
#' \describe{ \code{BFBayesFactorList} objects inherit from lists, and contain a
#' single slot:
#' 
#' \item{version}{character string giving the version and revision number of the
#' package that the model was created in}
#' 
#' Each element of the list contains a single
#' \code{"\link[=BFBayesFactor-class]{BFBayesFactor}"} object. Each element of
#' the list must have the same numerators, in the same order, as all the others.
#' The list object is displayed as a matrix of Bayes factors. }
#' @name BFBayesFactorList-class
#' @export
#' @examples
#' ## Compute some Bayes factors to demonstrate Bayes factor lists
#' data(puzzles)
#' bfs <- anovaBF(RT ~ shape*color + ID, data = puzzles, whichRandom = "ID", progress=FALSE)
#' 
#' ## Create a matrix of Bayes factors
#' bfList <- bfs / bfs
#' bfList
#' 
#' ## Use indexing to select parts of the 'matrix'
#' bfList[1,]
#' bfList[,1]
#' 
#' ## We can use the t (transpose) function as well, to get back a BFBayesFactor
#' t(bfList[2,])
#' 
#' ## Or transpose the whole matrix
#' t(bfList)
setClass("BFBayesFactorList", contains = "list", representation(version="character"))

#' @name BFBayesFactorTop-class
#' @rdname BFBayesFactor-class
setClass("BFBayesFactorTop", contains = "BFBayesFactor")

setOldClass("mcmc")
setClass("BFmcmc", contains = "mcmc", representation(model="BFmodel",data = "data.frame"))
