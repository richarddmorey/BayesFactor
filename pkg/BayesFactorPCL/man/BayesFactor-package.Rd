\name{BayesFactor-package}
\alias{BayesFactor-package}
\alias{BayesFactor}
\docType{package}
\title{
Functions to compute Bayes factor hypothesis tests for common research designs and hypotheses.
}
\description{
This package contains function to compute Bayes factors for a number of research designs and hypotheses, including
t tests, ANOVA, and linear regression. It also contains an experimental user interface.
}
\details{
\tabular{ll}{
Package: \tab BayesFactor\cr
Type: \tab Package\cr
Version: \tab 0.8.7\cr
Date: \tab 2012-11-22\cr
License: \tab GPL 2.0\cr
LazyLoad: \tab yes\cr
}
The following methods are currently implemented, with more to follow:

Linear regression: \code{\link{linearReg.Quad}},\code{\link{linearReg.Gibbs}} \code{\link{allLinearReg}}

t test: \code{\link{ttest.Quad}}, \code{\link{ttest.Gibbs}}

One-way ANOVA: \code{\link{oneWayAOV.Quad}}, \code{\link{oneWayAOV.Gibbs}}

General ANOVA: \code{\link{nWayAOV.MC}}, \code{\link{nWayAOV.Gibbs}}, \code{\link{allNways}}

Bayes factors for ANOVA designs can be analyzed using the experimental interface, \code{\link{aovGUI}}. 

}
\author{
Richard D. Morey and Jeffrey N. Rouder

Maintainer: Richard D. Morey <richarddmorey@gmail.com>
}
\references{
Liang, F. and Paulo, R. and Molina, G. and Clyde, M. A. and Berger, J. O. (2008). Mixtures of g-priors for Bayesian Variable Selection.
Journal of the American Statistical Association, 103, pp. 410-423

Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., \& Iverson, G. (2009). Bayesian t-tests
for accepting and rejecting the null hypothesis. Psychonomic Bulletin & Review, 16, 752-760

Rouder, J. N., Morey, R. D., Speckman, P. L., Province, J. M., (2012) Default Bayes Factors for ANOVA Designs. Journal of Mathematical Psychology.  56.  p. 356-374.

Perception and Cognition Lab (University of Missouri): Bayes factor calculators. \url{http://pcl.missouri.edu/bayesfactor}
}
\keyword{htest}
\seealso{
\code{\link[BAS:BAS-package]{BAS}} 
}
\examples{
## See specific functions for examples
}
