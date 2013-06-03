\name{generalTestBF}
\alias{generalTestBF}
\title{Function to compute Bayes factors for general designs}
\usage{
  generalTestBF(formula, data, whichRandom = NULL,
    denominator = NULL, whichModels = "withmain",
    neverExclude = NULL, iterations = 10000,
    progress = options()$BFprogress,
    rscaleFixed = "medium", rscaleRandom = "nuisance",
    rscaleCont = "medium", multicore = FALSE,
    method = "auto", noSample = FALSE)
}
\arguments{
  \item{formula}{a formula containing the full model for
  the analysis (see Examples)}

  \item{data}{a data frame containing data for all factors
  in the formula}

  \item{whichRandom}{a character vector specifying which
  factors are random}

  \item{denominator}{a (one-sided) formula containing the
  (optional) denominator for the analysis.  If present,
  will be added to each model. (see Note and Examples)}

  \item{whichModels}{which set of models to compare; see
  Details}

  \item{neverExclude}{a character vector containing a
  regular expression (see help for \link{regex} for
  details) that indicates which terms to always keep in the
  analysis}

  \item{iterations}{How many Monte Carlo simulations to
  generate, if relevant}

  \item{progress}{if \code{TRUE}, show progress with a text
  progress bar}

  \item{rscaleFixed}{prior scale for standardized, reduced
  fixed effects. A number of preset values can be given as
  strings; see Details.}

  \item{rscaleRandom}{prior scale for standardized random
  effects}

  \item{rscaleCont}{prior scale for standardized slopes}

  \item{multicore}{if \code{TRUE} use multiple cores
  through the \code{doMC} package. Unavailable on Windows.}

  \item{method}{approximation method, if needed. See
  \code{\link{nWayAOV}} for details.}

  \item{noSample}{if \code{TRUE}, do not sample, instead
  returning NA.}
}
\value{
  An object of class \code{BFBayesFactor}, containing the
  computed model comparisons
}
\description{
  This function computes Bayes factors corresponding to
  restrictions on a full model.
}
\details{
  See the help for \code{\link{anovaBF}} and
  \code{\link{anovaBF}} or details.

  Models, priors, and methods of computation are provided
  in Rouder et al. (2012) and Liang et al (2008).
}
\note{
  The function \code{generalTestBF} can compute Bayes
  factors for all restrictions of a full model against the
  null hypothesis that all effects are 0. The total number
  of tests computed -- if all tests are requested -- will
  be \eqn{2^K-1}{2^K - 1} for \eqn{K} factors or
  covariates.  This number increases very quickly with the
  number of tested predictors. An option is included to
  prevent testing too many models:
  \code{options('BFMaxModels')}, which defaults to 50,000,
  is the maximum number of models that will be analyzed at
  once. This can be increased by increased using
  \code{\link{options}}.

  It is possible to reduce the number of models tested by
  only testing the most complex model and every restriction
  that can be formed by removing one factor or interaction
  using the \code{whichModels} argument. See the help for
  \code{\link{anovaBF}} for details.

  To test the model against a specific \code{denominator},
  a (one-sided) formula can be specified which will be
  added to each model. These effects do not need to be
  specified in the \code{formula} argument. Most
  importantly, selection of models to be tested (which is
  controlled by arguments \code{whichModels} and
  \code{neverExclude}) is only performed on \code{formula}.
  Random effects can be specified for both \code{formula}
  and \code{denominator}. See Examples.
}
\examples{
## Puzzles example: see ?puzzles and ?anovaBF
data(puzzles)
## neverExclude argument makes sure that participant factor ID
## is in all models
result = generalTestBF(RT ~ shape*color + ID, data = puzzles, whichRandom = "ID",
neverExclude="ID", progress=FALSE)
result

## A more reasonable model could include random slopes in addition
## to the random intercept for ID
result2 = generalTestBF(RT ~ shape*color, data = puzzles, whichRandom = "ID",
denominator = ~ ID+ID:shape+ID:color, progress=FALSE)
result2
}
\author{
  Richard D. Morey (\email{richarddmorey@gmail.com}).
  Denominator argument added by Henrik Singmann.
}
\references{
  Rouder, J. N., Morey, R. D., Speckman, P. L., Province,
  J. M., (2012) Default Bayes Factors for ANOVA Designs.
  Journal of Mathematical Psychology.  56.  p. 356-374.

  Liang, F. and Paulo, R. and Molina, G. and Clyde, M. A.
  and Berger, J. O. (2008). Mixtures of g-priors for
  Bayesian Variable Selection. Journal of the American
  Statistical Association, 103, pp.  410-423
}
\seealso{
  \code{\link{lmBF}}, for testing specific models, and
  \code{\link{regressionBF}} and \code{anovaBF} for other
  functions for testing multiple models simultaneously.
}
\keyword{htest}

