\name{BayesFactorPCL-package}
\alias{BayesFactorPCL-package}
\alias{BayesFactorPCL}
\docType{package}
\title{
Functions to compute Bayes factor hypothesis tests for common research designs and hypotheses.
}
\description{
This package contains function to compute Bayes factors for a number of research designs and hypotheses, including
$t$ tests, ANOVA, and linear regression. 
}
\details{
\tabular{ll}{
Package: \tab BayesFactorPCL\cr
Type: \tab Package\cr
Version: \tab 0.6\cr
Date: \tab 2010-11-03\cr
License: \tab GPL 2.0\cr
LazyLoad: \tab yes\cr
}
The following functions are currently implemented, with more to follow:

ttest.Quad(): Computes a Bayes factor for one- and two-sample designs

oneWayAOV.Quad(): Computes a Bayes factor for one-way ANOVA-type designs
}
\author{
Richard D. Morey and Jeffrey N. Rouder

Maintainer: Richard D. Morey <richarddmorey@gmail.com>
}
\references{
Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., \& Iverson, G. (2009). Bayesian t-tests
for accepting and rejecting the null hypothesis. Psychonomic Bulletin & Review, 16, 752-760

Morey, R. D., Rouder, J. N., Pratte, M. S., \& Speckman, P. L. (submitted). Conditional Marginal Density 
Estimation to Improve Bayes Factor Computation.

Perception and Cognition Lab (University of Missouri): Bayes factor calculators. \url{http://pcl.missouri.edu/bayesfactor}
}
\keyword{htest}
\seealso{
\code{\link[BAS:BAS-package]{BAS}} ~~
}
\examples{
## See specific functions for examples
}
