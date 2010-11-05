\name{ttest.Quad}
\alias{ttest.Quad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to compute the Bayes factors for one- or two-sample designs
}
\description{
This function computes the Bayes factor corresponding to a one-sample or two-sample $t$~test, using Gaussian quadrature integration. 
See details about priors in the details.
}
\usage{
ttest.Quad(t, n1, n2 = 0, rscale = 1, prior.cauchy = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
	The $t$ value obtained from the corresponding $t$~test
}
  \item{n1}{
Sample size for first group
}
  \item{n2}{
Sample size for second group. If a one-sample test is desired, \code{n1} should be 0.
}
  \item{rscale}{
The prior scale (see Details below)
}
  \item{prior.cauchy}{
Use a Cauchy prior for the prior on effect size? If false, the prior on effect size is Normal (see Details below) 
}
}
\details{
This Bayes factor tests the null hypothesis that the mean of a normal population is 0. 
Specifically, the Bayes factor compares two hypotheses: that the standardized effect size is 0, or that the 
standardized effect size is not 0. For one-sample tests, the standardized effect size is $\mu/\sigma$; 
for two sample tests, the standardized effect size is $(\mu_2-\mu_1)/\sigma$.

A noninformative Jeffreys prior is placed on the variance of the normal population, while a Cauchy 
prior is placed on the standardized effect size. The \code{rscale} argument controls the scale of the prior 
distribution, with \code{rscale=1} yielding a standard Cauchy prior. If \code{prior.Cauchy=FALSE}, then a Normal prior is used
instead of a Cauchy prior. See the references below for more details.
}
\value{
A scalar giving the Bayes factor in favor of the null hypothesis that the effect size is 0.
}
\references{

\item{}{Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., \& Iverson, G. (2009). Bayesian t–tests
for accepting and rejecting the null hypothesis. Psychonomic Bulletin & Review, 16, 752-760}

\item{}{Perception and Cognition Lab (University of Missouri): Bayes factor calculators. \link{http://pcl.missouri.edu/bayesfactor}}
}
\author{
Jeffrey N. Rouder and Richard D. Morey
}
\note{
For samples from the marginal posterior distributions of the parameters, see \code{\link{ttest.Gibbs}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{t.test}}, \code{\link{ttest.Gibbs}}, \code{\link{oneWayAOV.Quad}},\code{\link{oneWayAOV.Gibbs}}
}
\examples{
## Classical example, taken from t.test() example
## Student's sleep data
plot(extra ~ group, data = sleep)

## t.test() gives a t value of -1.8608
t.test(extra ~ group, data = sleep, var.equal=TRUE)

## Gives a Bayes factor of about 1.02
## in favor of the null hypothesis
ttest.Quad(t = -1.8608, n1 = 10, n2 = 10, rscale=sqrt(2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
