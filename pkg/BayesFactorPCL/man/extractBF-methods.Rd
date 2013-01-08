\docType{methods}
\name{extractBF}
\alias{extractBF}
\alias{extractBF,BFBayesFactor-method}
\title{Extract the Bayes factor from an object}
\arguments{
  \item{x}{object from which to extract the Bayes factor}

  \item{...}{arguments passed to and from related methods}
}
\value{
  Returns an object containing Bayes factors extracted from
  the object
}
\description{
  Extract the Bayes factor from an object
}
\examples{
## Sample from the posteriors for two models
data(puzzles)

bf = lmBF(RT ~ shape*color + ID, data = puzzles, whichRandom="ID", progress=FALSE)

extractBF(bf)
}

