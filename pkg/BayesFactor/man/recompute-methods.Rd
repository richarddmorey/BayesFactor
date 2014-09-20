\docType{methods}
\name{recompute}
\alias{recompute}
\alias{recompute,BFBayesFactor-method}
\alias{recompute,BFmcmc-method}
\title{Take an object and redo the computation (useful for sampling)}
\usage{
recompute(x, progress = options()$BFprogress, ...)

\S4method{recompute}{BFBayesFactor}(x, progress = options()$BFprogress, ...)

\S4method{recompute}{BFmcmc}(x, progress, ...)
}
\arguments{
  \item{x}{object to recompute}

  \item{progress}{report progress of the computation?}

  \item{...}{arguments passed to and from related methods}
}
\value{
Returns an object of the same type, after repeating the
sampling (perhaps with more iterations)
}
\description{
Take an object and redo the computation (useful for
sampling)
}
\examples{
## Sample from the posteriors for two models
data(puzzles)

## Main effects model; result is a BFmcmc object, inheriting
## mcmc from the coda package
bf = lmBF(RT ~ shape + color + ID, data = puzzles, whichRandom = "ID",
   progress = FALSE)

## recompute Bayes factor object
recompute(bf, iterations = 1000, progress = FALSE)

## Sample from posterior distribution of model above, and recompute:
chains = posterior(bf, iterations = 1000, progress = FALSE)
newChains = recompute(chains, iterations = 1000, progress=FALSE)
}

