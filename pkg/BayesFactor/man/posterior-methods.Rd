\docType{methods}
\name{posterior}
\alias{posterior}
\alias{posterior,BFBayesFactor,missing,missing,missing-method}
\alias{posterior,BFBayesFactor,missing,missing,numeric-method}
\alias{posterior,BFBayesFactor,numeric,missing,numeric-method}
\alias{posterior,BFindepSample,missing,data.frame,numeric-method}
\alias{posterior,BFlinearModel,missing,data.frame,missing-method}
\alias{posterior,BFlinearModel,missing,data.frame,numeric-method}
\alias{posterior,BFoneSample,missing,data.frame,numeric-method}
\title{Sample from the posterior distribution of one of several models.}
\arguments{
  \item{model}{or set of models from which to sample}

  \item{index}{the index within the set of models giving
  the desired model}

  \item{data}{the data to be conditioned on}

  \item{iterations}{the number of iterations to sample}

  \item{...}{arguments passed to and from related methods}
}
\value{
  Returns an object containing samples from the posterior
  distribution of the specified model
}
\description{
  This function samples from the posterior distribution of
  a \code{BFmodel}, which can be obtained from a
  \code{BFBayesFactor} object. If there is more than one
  numerator in the \code{BFBayesFactor} object, the
  \code{index} argument can be passed to select one
  numerator.
}
\details{
  The data argument is used internally, and will typically
  not be needed by end-users.

  Note that if there are fixed effects in the model, the
  reduced parameterzation used internally (see help for
  \code{\link{anovaBF}}) is unreduced. For a factor with
  two levels, the chain will contain two effect estimates
  that sum to 0.

  Two useful arguments that can be passed to related
  methods are \code{thin} and \code{columnFilter},
  currently implemented for methods using \code{nWayAOV}
  (models with more than one categorical covariate, or a
  mix of categorical and continuous covariates).
  \code{thin}, an integer, will keep only every \code{thin}
  iterations. The default is \code{thin=1}, which keeps all
  iterations. Argument \code{columnFilter} is either
  \code{NULL} (for no filtering) or a character vector of
  extended regular expressions (see \link{regex} help for
  details). Any column from an effect that matches one of
  the filters will not be saved.
}
\examples{
## Sample from the posteriors for two models
data(sleep)

bf = lmBF(extra ~ group + ID, data = sleep, whichRandom="ID", progress=FALSE)

## sample from the posterior of the numerator model
## data argument not needed - it is included in the Bayes factor object
chains = posterior(bf, iterations = 1000, progress = FALSE)

plot(chains)

## demonstrate column filtering by filtering out participant effects
data(puzzles)
bf = lmBF(RT ~ shape + color + shape:color + ID, data=puzzles)
chains = posterior(bf, iterations = 1000, progress = FALSE, columnFilter="^ID$")
colnames(chains) # Contains no participant effects
}

