\name{oneWayAOV.Fstat}
\alias{oneWayAOV.Fstat}
\title{Use F statistic to compute Bayes factor for balanced one-way designs}
\usage{
  oneWayAOV.Fstat(F, N, J, rscale = "medium")
}
\arguments{
  \item{F}{F statistic from classical ANOVA}

  \item{N}{number of observations per cell or group}

  \item{J}{number of cells or groups}

  \item{rscale}{numeric prior scale}
}
\value{
  a vector of length 2 containing the computed log(e) Bayes
  factor (against the intercept-only null), along with a
  proportional error estimate on the Bayes factor.
}
\description{
  Using the classical F test statistic for a balanced
  one-way design, this function computes the corresponding
  Bayes factor test.
}
\details{
  For F statistics computed from balanced one-way designs,
  this function can be used to compute the Bayes factor
  testing the model that all group means are not equal to
  the grand mean, versus the null model that all group
  means are equal. It can be used when you don't have
  access to the full data set for analysis by
  \code{\link{lmBF}}, but you do have the test statistic.

  For details about the model, see the help for
  \code{\link{anovaBF}}, and the references therein.

  The Bayes factor is computed via Gaussian quadrature.
}
\note{
  \code{oneWayAOV.Fstat} should only be used with F values
  obtained from balanced designs.
}
\examples{
## Example data "InsectSprays" - see ?InsectSprays
require(stats); require(graphics)
boxplot(count ~ spray, data = InsectSprays, xlab = "Type of spray", ylab = "Insect count", main = "InsectSprays data", varwidth = TRUE, col = "lightgray")

## Classical analysis (with transformation)
classical <- aov(sqrt(count) ~ spray, data = InsectSprays)
plot(classical)
summary(classical)

## Bayes factor (a very large number)
Fvalue <- anova(classical)$"F value"[1]
result <- oneWayAOV.Fstat(Fvalue, N=12, J=6)
exp(result[['bf']])
}
\author{
  Richard D. Morey (\email{richarddmorey@gmail.com})
}
\references{
  Morey, R. D., Rouder, J. N., Pratte, M. S., \& Speckman,
  P. L.  (2011). Using MCMC chain outputs to efficiently
  estimate Bayes factors.  Journal of Mathematical
  Psychology, 55, 368-378
}
\seealso{
  \code{\link{integrate}}, \code{\link{aov}}; see
  \code{\link{lmBF}} for the intended interface to this
  function, using the full data set.
}
\keyword{htest}

