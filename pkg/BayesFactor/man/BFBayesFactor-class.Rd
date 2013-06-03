\name{BFBayesFactor-class}
\alias{/,BFBayesFactor,BFBayesFactor-method}
\alias{/,numeric,BFBayesFactor-method}
\alias{[,BFBayesFactor,index,missing,missing-method}
\alias{[,BFBayesFactorTop,index,missing,missing-method}
\alias{BFBayesFactor-class}
\alias{BFBayesFactorTop-class}
\alias{is.na,BFBayesFactor-method}
\alias{t,BFBayesFactor-method}
\alias{which.max,BFBayesFactor-method}
\alias{which.min,BFBayesFactor-method}
\title{General S4 class for representing multiple Bayes factor model comparisons, all against the same model}
\description{
  The \code{BFBayesFactor} class is a general S4 class for
  representing models model comparison via Bayes factor.
}
\details{
  \code{BFBayesFactor} objects can be inverted by taking
  the reciprocal and can be divided by one another,
  provided both objects have the same denominator. In
  addition, the \code{t} (transpose) method can be used to
  invert Bayes factor objects.  \describe{ The
  \code{BFBayesFactor} class has the following slots
  defined: \item{numerator}{a list of models all inheriting
  \code{BFmodel}, each providing a single denominator}
  \item{denominator}{a single \code{BFmodel} object serving
  as the denominator for all model comparisons}
  \item{bayesFactor}{a data frame containing information
  about the comparison between each numerator and the
  denominator} \item{data}{a data frame containing the data
  used for the comparison} \item{version}{character string
  giving the version and revision number of the package
  that the model was created in} }
}
\examples{
## Compute some Bayes factors to demonstrate division and indexing
data(puzzles)
bfs <- anovaBF(RT ~ shape*color + ID, data = puzzles, whichRandom = "ID", progress=FALSE)

## First and second models can be separated; they remain BFBayesFactor objects
b1 = bfs[1]
b2 = bfs[2]
b1

## We can invert them, or divide them to obtain new model comparisons
1/b1
b1 / b2

## Use transpose to create a BFBayesFactorList
t(bfs)
}

