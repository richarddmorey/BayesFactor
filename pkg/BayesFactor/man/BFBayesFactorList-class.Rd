\docType{class}
\name{BFBayesFactorList-class}
\alias{/,numeric,BFBayesFactorList-method}
\alias{BFBayesFactorList-class}
\alias{[,BFBayesFactorList,index,index,missing-method}
\alias{[,BFBayesFactorList,index,missing,missing-method}
\alias{[,BFBayesFactorList,missing,index,missing-method}
\alias{t,BFBayesFactorList-method}
\title{General S4 class for representing a collection of Bayes factor model
comprisons, each against a different denominator}
\usage{
\S4method{t}{BFBayesFactorList}(x)

\S4method{/}{numeric,BFBayesFactorList}(e1, e2)

\S4method{[}{BFBayesFactorList,index,index,missing}(x, i, j, ..., drop)

\S4method{[}{BFBayesFactorList,index,missing,missing}(x, i, j, ..., drop)

\S4method{[}{BFBayesFactorList,missing,index,missing}(x, i, j, ..., drop)
}
\arguments{
  \item{x}{a matrix or data frame, typically.}

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
The \code{BFBayesFactorList} class is a general S4 class
for representing models model comparison via Bayes factor.
See the examples for demonstrations of BFBayesFactorList
methods.
}
\details{
\describe{ \code{BFBayesFactorList} objects inherit from
lists, and contain a single slot:

\item{version}{character string giving the version and
revision number of the package that the model was created
in}

Each element of the list contains a single
\code{"\link[=BFBayesFactor-class]{BFBayesFactor}"} object.
Each element of the list must have the same numerators, in
the same order, as all the others. The list object is
displayed as a matrix of Bayes factors. }
}
\examples{
## Compute some Bayes factors to demonstrate Bayes factor lists
data(puzzles)
bfs <- anovaBF(RT ~ shape*color + ID, data = puzzles, whichRandom = "ID", progress=FALSE)

## Create a matrix of Bayes factors
bfList <- bfs / bfs
bfList

## Use indexing to select parts of the 'matrix'
bfList[1,]
bfList[,1]

## We can use the t (transpose) function as well, to get back a BFBayesFactor
t(bfList[2,])

## Or transpose the whole matrix
t(bfList)
}

