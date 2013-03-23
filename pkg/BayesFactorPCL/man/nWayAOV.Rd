\name{nWayAOV}
\alias{nWayAOV}
\title{Use ANOVA design matrix to compute Bayes factors or sample posterior}
\usage{
  nWayAOV(y, X, struc = NULL, gMap = NULL, rscale,
    iterations = 10000, progress = TRUE, gibi = NULL,
    gibbs = FALSE, ignoreCols = NULL, thin = 1,
    method = "auto", continuous = FALSE)
}
\arguments{
  \item{y}{vector of observations}

  \item{X}{design matrix whose number of rows match
  \code{length(y)}.}

  \item{struc}{vector grouping the columns of \code{X} (see
  Details).}

  \item{gMap}{alternative way of grouping the columns of
  \code{X}}

  \item{rscale}{a vector of prior scale(s) of appropriate
  length (see Details).}

  \item{iterations}{Number of Monte Carlo samples used to
  estimate Bayes factor or posterior}

  \item{progress}{if \code{TRUE}, show progress with a text
  progress bar}

  \item{gibi}{interface for a future graphical user
  interface (not intended for use by end users)}

  \item{gibbs}{if \code{TRUE}, return samples from the
  posterior instead of a Bayes factor}

  \item{ignoreCols}{if \code{NULL} and \code{gibbs=TRUE},
  all parameter estimates are returned in the MCMC object.
  If not \code{NULL}, a vector of length P-1 (where P is
  number of columns in the design matrix) giving which
  effect estimates to ignore in output}

  \item{thin}{MCMC chain to every \code{thin} iterations.
  Default of 1 means no thinning. Only used if
  \code{gibbs=TRUE}}

  \item{method}{the integration method (only valid if
  \code{gibbs=TRUE}); one of "simple", "importance",
  "laplace", or "auto"}

  \item{continuous}{either FALSE is no continuous
  covariates are included, or a logical vector of length
  equal to number of columns of X indicating which columns
  of the design matrix represent continuous covariates}
}
\value{
  If \code{posterior} is \code{FALSE}, a vector of length 2
  containing the computed log(e) Bayes factor (against the
  intercept-only null), along with a proportional error
  estimate on the Bayes factor. Otherwise, an object of
  class \code{mcmc}, containing MCMC samples from the
  posterior is returned.
}
\description{
  Computes a single Bayes factor, or samples from the
  posterior, for an ANOVA model defined by a design matrix
}
\details{
  This function is not meant to be called by end-users,
  although technically-minded users can call this function
  for flexibility beyond what the other functions in this
  package provide. See \code{\link{lmBF}} for a
  user-friendly front-end to this function. Details about
  the priors can be found in the help for
  \code{\link{anovaBF}} and the references therein.

  Arguments \code{struc} and \code{gMap} provide a way of
  grouping columns of the design matrix as a factor; the
  effects in each group will share a common \eqn{g}
  parameter. Only one of these arguments is needed; if both
  are given, \code{gMap} takes precedence.

  \code{gMap} should be a vector of the same length as the
  number of nonconstant rows in \code{X}. It will contain
  all integers from 0 to \eqn{N_g-1}{Ng-1}, where
  \eqn{N_g}{Ng} is the total number of \eqn{g} parameters.
  Each element of \code{gMap} specifies the group to which
  that column belongs.

  If all columns belonging to a group are adjacent,
  \code{struc} can instead be used to compactly represent
  the groupings. \code{struc} is a vector of length
  \eqn{N_g}{Ng}. Each element specifies the number columns
  in the group. \code{gMap} is thus the
  \code{\link{inverse.rle}} of \code{struc}, minus 1.

  The vector \code{rscale} should be of length
  \eqn{N_g}{Ng}, and contain the prior scales of the
  standardized effects. See Rouder et al. (2012) for more
  details and the help for \code{\link{anovaBF}} for some
  typical values.

  The method used to estimate the Bayes factor depends on
  the \code{method} argument. "simple" is most accurate for
  small to moderate sample sizes, and uses the Monte Carlo
  sampling method described in Rouder et al. (2012).
  "importance" uses an importance sampling algorithm with
  an importance distribution that is multivariate normal on
  log(g). "laplace" does not sample, but uses a Laplace
  approximation to the integral. It is expected to be more
  accurate for large sample sizes, where MC sampling is
  slow. If \code{method="auto"}, then an initial run with
  both samplers is done, and the sampling method that
  yields the least-variable samples is chosen. The number
  of initial test iterations is determined by
  \code{options(BFpretestIterations)}.

  If posterior samples are requested, the posterior is
  sampled with a Gibbs sampler.
}
\examples{
## Classical example, taken from t.test() example
## Student's sleep data
data(sleep)
plot(extra ~ group, data = sleep)

## traditional ANOVA gives a p value of 0.00283
summary(aov(extra ~ group + Error(ID/group), data = sleep))

## Build design matrix
group.column <- rep(1/c(-sqrt(2),sqrt(2)),each=10)
subject.matrix <- model.matrix(~sleep$ID - 1,data=sleep$ID)
## Note that we include no constant column
X <- cbind(group.column, subject.matrix)

## (log) Bayes factor of full model against grand-mean only model
bf.full <- nWayAOV(y = sleep$extra, X = X, struc = c(1,10), rscale=c(.5,1))
exp(bf.full['bf'])

## Compare with lmBF result (should be about the same, give or take 1\%)
bf.full2 <- lmBF(extra ~ group + ID, data = sleep, whichRandom = "ID")
bf.full2
}
\author{
  Richard D. Morey (\email{richarddmorey@gmail.com}),
  Jeffery N.  Rouder (\email{rouderj@missouri.edu})
}
\references{
  Rouder, J. N., Morey, R. D., Speckman, P. L., Province,
  J. M., (2012) Default Bayes Factors for ANOVA Designs.
  Journal of Mathematical Psychology.  56.  p. 356-374.
}
\seealso{
  See \code{\link{lmBF}} for the user-friendly front end to
  this function; see \code{\link{regressionBF}} and
  \code{anovaBF} for testing many regression or ANOVA
  models simultaneously.
}
\keyword{htest}

