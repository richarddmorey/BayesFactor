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
  \item{x}{a BFBayesFactorList object}

  \item{e1}{Numerator of the ratio}

  \item{e2}{Denominator of the ratio}

  \item{i}{indices specifying rows to extract}

  \item{j}{indices specifying columns to extract}

  \item{drop}{unused}

  \item{...}{further arguments passed to related methods}
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

