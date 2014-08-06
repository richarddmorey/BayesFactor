\name{meta.ttest.tstat}
\alias{meta.ttest.tstat}
\title{Use t statistic to compute Bayes factor for one- and two- sample designs}
\usage{
meta.ttest.tstat(t, n1, n2 = NULL, nullInterval = NULL, rscale = "medium")
}
\arguments{
  \item{t}{vector of classical t statistics}

  \item{n1}{vector of sizes of first group (or only group,
  for one-sample tests)}

  \item{n2}{vector of size of second group, for
  independent-groups tests}

  \item{nullInterval}{optional vector of length 2
  containing lower and upper bounds of an interval
  hypothesis to test, in standardized units}

  \item{rscale}{numeric prior scale}
}
\value{
If \code{nullInterval} is defined, then two Bayes factors
will be computed: The Bayes factor for the interval against
the null hypothesis that the standardized effect is 0, and
the corresponding Bayes factor for the compliment of the
interval. For each Bayes factor, a vector of length 2
containing the computed log(e) Bayes factor (against the
point null), along with a proportional error estimate on
the Bayes factor is returned.
}
\description{
Using the classical t test statistic for multiple one- or
two-sample designed-experiments, this function computes the
corresponding meta-analytic Bayes factor test.
}
\details{
This function can be used to compute a meta-analytic Bayes
factor corresponding to a one-sample, a paired-sample, or
an independent-groups t test, using the classical t
statistic.

For details about the model, see the references below, and
the help for \code{\link{ttestBF}}, and the references
therein.

The Bayes factor is computed via Gaussian quadrature.
}
\examples{
## Bem's (2010) data (see Rouder & Morey, 2011)
t=c(-.15,2.39,2.42,2.43)
N=c(100,150,97,99)

## Using rscale=1 to be consistent with Rouder & Morey (2011)
exp(meta.ttest.tstat(t,n1=N,rscale=1)$bf)
}
\author{
Richard D. Morey (\email{richarddmorey@gmail.com})
}
\references{
Morey, R. D. & Rouder, J. N. (2011). Bayes Factor
Approaches for Testing Interval Null Hypotheses.
Psychological Methods, 16, 406-419

Rouder, J. N. & Morey, R. D. (2011). A Bayes Factor
Meta-Analysis of Bem's ESP Claim.  Psychonomic Bulletin &
Review, 18, 682-689
}
\seealso{
\code{\link{integrate}}, \code{\link{t.test}}; see
\code{\link{ttest.tstat}} for the single-study version of
this function
}
\keyword{htest}

