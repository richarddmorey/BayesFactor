\name{BayesFactor-package}
\alias{BayesFactor-package}
\alias{BayesFactor}
\docType{package}
\title{
Functions to compute Bayes factor hypothesis tests for common research designs and hypotheses.
}
\description{
This package contains function to compute Bayes factors for a number of research designs and hypotheses, including
t tests, ANOVA, and linear regression. 
}
\details{
\tabular{ll}{
Package: \tab BayesFactor\cr
Type: \tab Package\cr
Version: \tab 0.8.3\cr
Date: \tab 2012-10-11\cr
License: \tab GPL 2.0\cr
LazyLoad: \tab yes\cr
}
The following methods are currently implemented, with more to follow:

Linear regression: \code{\link{linearReg.Quad}}

t test: \code{\link{ttest.Quad}}, \code{\link{ttest.Gibbs}}

One-way ANOVA: \code{\link{oneWayAOV.Quad}}, \code{\link{oneWayAOV.Gibbs}}

General ANOVA: \code{\link{nWayAOV.MC}}, \code{\link{nWayAOV.Gibbs}}, \code{\link{allNways}}

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
\code{\link[BAS:BAS-package]{BAS}} 
}
\examples{
## See specific functions for examples
}
