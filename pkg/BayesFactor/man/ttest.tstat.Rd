\name{ttest.tstat}
\alias{ttest.tstat}
\title{Use t statistic to compute Bayes factor for one- and two- sample designs}
\usage{
  ttest.tstat(t, n1, n2 = 0, nullInterval = NULL,
    rscale = "medium")
}
\arguments{
  \item{t}{classical t statistic}

  \item{n1}{size of first group (or only group, for
  one-sample tests)}

  \item{n2}{size of second group, for independent-groups
  tests}

  \item{nullInterval}{optional vector of length 2
  containing lower and upper bounds of an interval
  hypothesis to test, in standardized units}

  \item{rscale}{numeric prior scale}
}
\value{
  If \code{nullInterval} is defined, then two Bayes factors
  will be computed: The Bayes factor for the interval
  against the null hypothesis that the standardized effect
  is 0, and the corresponding Bayes factor for the
  compliment of the interval. For each Bayes factor, a
  vector of length 2 containing the computed log(e) Bayes
  factor (against the point null), along with a
  proportional error estimate on the Bayes factor is
  returned.
}
\description{
  Using the classical t test statistic for a one- or
  two-sample design, this function computes the
  corresponding Bayes factor test.
}
\details{
  This function can be used to compute the Bayes factor
  corresponding to a one-sample, a paired-sample, or an
  independent-groups t test, using the classical t
  statistic.  It can be used when you don't have access to
  the full data set for analysis by \code{\link{ttestBF}},
  but you do have the test statistic.

  For details about the model, see the help for
  \code{\link{ttestBF}}, and the references therein.

  The Bayes factor is computed via Gaussian quadrature.
}
\examples{
## Classical example: Student's sleep data
data(sleep)
plot(extra ~ group, data = sleep)

## t.test() gives a t value of -4.0621
t.test(extra ~ group, data = sleep, paired=TRUE)
## Gives a Bayes factor of about 15
## in favor of the alternative hypothesis
result <- ttest.tstat(t = -4.0621, n1 = 10)
exp(result[['bf']])
}
\author{
  Richard D. Morey (\email{richarddmorey@gmail.com}) and
  Jeffrey N.  Rouder (\email{rouderj@missouri.edu})
}
\references{
  Morey, R. D. & Rouder, J. N. (2011). Bayes Factor
  Approaches for Testing Interval Null Hypotheses.
  Psychological Methods, 16, 406-419

  Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., &
  Iverson, G.  (2009). Bayesian t-tests for accepting and
  rejecting the null hypothesis.  Psychonomic Bulletin &
  Review, 16, 752-760
}
\seealso{
  \code{\link{integrate}}, \code{\link{t.test}}; see
  \code{\link{ttestBF}} for the intended interface to this
  function, using the full data set.
}
\keyword{htest}

