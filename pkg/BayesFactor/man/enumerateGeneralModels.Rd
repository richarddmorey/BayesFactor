\name{enumerateGeneralModels}
\alias{enumerateGeneralModels}
\title{Function for generation of nested linear models}
\usage{
enumerateGeneralModels(fmla, whichModels, neverExclude = NULL,
  includeBottom = TRUE)
}
\arguments{
  \item{fmla}{formula for the "full" model}

  \item{whichModels}{which subsets of models to generate}

  \item{neverExclude}{a character vector of terms to never
  remove}

  \item{includeBottom}{Include the base model containing
  only \code{neverExclude} terms}
}
\description{
Generate lists of nested models, given a model formula
}
\details{
This is a backend function not intended for users. It is
exposed for third-party applications.
}
\keyword{internal}

