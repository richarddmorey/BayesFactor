\docType{methods}
\name{compare}
\alias{compare}
\alias{compare,BFcontingencyTable,BFcontingencyTable,data.frame-method}
\alias{compare,BFcontingencyTable,missing,data.frame-method}
\alias{compare,BFindepSample,missing,data.frame-method}
\alias{compare,BFlinearModel,BFlinearModel,data.frame-method}
\alias{compare,BFlinearModel,missing,data.frame-method}
\alias{compare,BFmcmc,BFmcmc,ANY-method}
\alias{compare,BFmcmc,missing,ANY-method}
\alias{compare,BFoneSample,missing,data.frame-method}
\title{Compare two models, with respect to some data}
\usage{
compare(numerator, denominator, data, ...)
}
\arguments{
  \item{numerator}{first model}

  \item{denominator}{second model (if omitted, compare to
  predefined null)}

  \item{data}{data for the comparison}

  \item{...}{arguments passed to and from related methods}
}
\value{
The compare function will return a model comparison object,
typically a Bayes factor
}
\description{
This method is used primarily in the backend, and will only
rarely be called by the end user. But see the examples
below for a demonstration.
}
\examples{
## Sample from the posteriors for two models
data(puzzles)

## Main effects model; result is a BFmcmc object, inheriting
## mcmc from the coda package
mod1 = lmBF(RT ~ shape + color + ID, data = puzzles, whichRandom = "ID",
   progress = FALSE, posterior = TRUE, iterations = 1000)

plot(mod1)

## Full model
mod2 = lmBF(RT ~ shape*color + ID, data = puzzles, whichRandom = "ID",
   progress = FALSE, posterior = TRUE, iterations = 1000)

## Each BFmcmc object contains the model used to generate it, so we
## can compare them (data is not needed, it is contained in the objects):

compare(mod1, mod2)
}

