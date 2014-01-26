\docType{class}
\name{BFBayesFactor-class}
\alias{/,BFBayesFactor,BFBayesFactor-method}
\alias{/,numeric,BFBayesFactor-method}
\alias{BFBayesFactor-class}
\alias{BFBayesFactorTop-class}
\alias{[,BFBayesFactor,index,missing,missing-method}
\alias{[,BFBayesFactorTop,index,missing,missing-method}
\alias{is.na,BFBayesFactor-method}
\alias{t,BFBayesFactor-method}
\alias{which.max,BFBayesFactor-method}
\alias{which.min,BFBayesFactor-method}
\title{General S4 class for representing multiple Bayes factor model comparisons, all against the same model}
\usage{
\S4method{/}{numeric,BFBayesFactor}(e1, e2)

\S4method{/}{BFBayesFactor,BFBayesFactor}(e1, e2)

\S4method{[}{BFBayesFactor,index,missing,missing}(x, i, j, ..., drop)

\S4method{t}{BFBayesFactor}(x)

\S4method{which.max}{BFBayesFactor}(x)

\S4method{which.min}{BFBayesFactor}(x)

\S4method{is.na}{BFBayesFactor}(x)

\S4method{[}{BFBayesFactorTop,index,missing,missing}(x, i, j, ..., drop)
}
\arguments{
  \item{e1}{Numerator of the ratio}

  \item{e2}{Denominator of the ratio}

  \item{j}{unused for BFBayesFactor objects}

  \item{...}{further arguments passed to related methods}

  \item{x}{numeric or complex vectors or objects which can
  be coerced to such, or other objects for which methods
  have been written.}

  \item{i}{ indices specifying elements to extract or
  replace.  Indices are \code{numeric} or \code{character}
  vectors or empty (missing) or \code{NULL}.  Numeric
  values are coerced to integer as by
  \code{\link{as.integer}} (and hence truncated towards
  zero).  Character vectors will be matched to the
  \code{\link{names}} of the object (or for
  matrices/arrays, the \code{\link{dimnames}}): see
  \sQuote{Character indices} below for further details.

  For \code{[}-indexing only: \code{i}, \code{j},
  \code{\dots} can be logical vectors, indicating
  elements/slices to select.  Such vectors are recycled if
  necessary to match the corresponding extent.  \code{i},
  \code{j}, \code{\dots} can also be negative integers,
  indicating elements/slices to leave out of the selection.

  When indexing arrays by \code{[} a single argument
  \code{i} can be a matrix with as many columns as there
  are dimensions of \code{x}; the result is then a vector
  with elements corresponding to the sets of indices in
  each row of \code{i}.

  An index value of \code{NULL} is treated as if it were
  \code{integer(0)}.  }

  \item{drop}{For matrices and arrays.  If \code{TRUE} the
  result is coerced to the lowest possible dimension (see
  the examples).  This only works for extracting elements,
  not for the replacement.  See \code{\link{drop}} for
  further details.  }
}
\description{
The \code{BFBayesFactor} class is a general S4 class for
representing models model comparison via Bayes factor.
}
\details{
\code{BFBayesFactor} objects can be inverted by taking the
reciprocal and can be divided by one another, provided both
objects have the same denominator. In addition, the
\code{t} (transpose) method can be used to invert Bayes
factor objects.  \describe{ The \code{BFBayesFactor} class
has the following slots defined: \item{numerator}{a list of
models all inheriting \code{BFmodel}, each providing a
single denominator} \item{denominator}{a single
\code{BFmodel} object serving as the denominator for all
model comparisons} \item{bayesFactor}{a data frame
containing information about the comparison between each
numerator and the denominator} \item{data}{a data frame
containing the data used for the comparison}
\item{version}{character string giving the version and
revision number of the package that the model was created
in} }
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

