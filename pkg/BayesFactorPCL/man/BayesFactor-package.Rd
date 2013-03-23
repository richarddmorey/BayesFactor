\docType{package}
\name{BayesFactor-package}
\alias{BayesFactor}
\alias{BayesFactor-package}
\title{Functions to compute Bayes factor hypothesis tests for common research designs
and hypotheses.}
\description{
  This package contains function to compute Bayes factors
  for a number of research designs and hypotheses,
  including t tests, ANOVA, and linear regression.
}
\details{
  \tabular{ll}{ Package: \tab BayesFactor\cr Type: \tab
  Package\cr Version: \tab 0.9.5\cr Date: \tab 2013-3-20\cr
  License: \tab GPL 2.0\cr LazyLoad: \tab yes\cr } The
  following methods are currently implemented, with more to
  follow:

  Linear regression: \code{\link{regressionBF}}
  \code{\link{lmBF}}, \code{\link{linearReg.R2stat}};

  t test: \code{\link{ttestBF}}, \code{\link{ttest.tstat}};

  ANOVA: \code{\link{anovaBF}}, \code{\link{lmBF}},
  \code{\link{oneWayAOV.Fstat}};

  Other useful functions: \code{\link{posterior}}, for
  sampling from posterior distributions;
  \code{\link{recompute}}, for re-estimating a Bayes factor
  or posterior distribution; \code{\link{compare}}, to
  compare two model posteriors; and
  \code{\link{plot.BFBayesFactor}}, for plotting Bayes
  factor objects.
}
\examples{
## See specific functions for examples.
}
\author{
  Richard D. Morey and Jeffrey N. Rouder

  Maintainer: Richard D. Morey <richarddmorey@gmail.com>
}
\references{
  Liang, F. and Paulo, R. and Molina, G. and Clyde, M. A.
  and Berger, J. O. (2008). Mixtures of g-priors for
  Bayesian Variable Selection.  Journal of the American
  Statistical Association, 103, pp. 410-423

  Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., \&
  Iverson, G.  (2009). Bayesian t-tests for accepting and
  rejecting the null hypothesis.  Psychonomic Bulletin &
  Review, 16, 752-760

  Rouder, J. N., Morey, R. D., Speckman, P. L., Province,
  J. M., (2012) Default Bayes Factors for ANOVA Designs.
  Journal of Mathematical Psychology.  56.  p. 356-374.

  Perception and Cognition Lab (University of Missouri):
  Bayes factor calculators.
  \url{http://pcl.missouri.edu/bayesfactor}
}
\seealso{
  \code{\link[BAS:BAS-package]{BAS}}
}
\keyword{htest}

