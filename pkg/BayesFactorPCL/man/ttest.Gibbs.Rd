\name{t tests}
\alias{ttest.Gibbs}
\alias{ttest.Quad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to for Bayesian analysis of one- and two sample designs 
}
\description{
\code{ttest.Gibbs} samples from the posterior distributions of parameters for 
one-sample designs, using Gibbs sampling; \code{ttest.Quad} computes Bayes factors for one- or two-sample designs.
}
\usage{
ttest.Quad(t, n1, n2 = 0, rscale = "medium", prior.cauchy = TRUE, 
    logbf=FALSE)

ttest.Gibbs(y, t, n, iterations=10000, rscale="medium", null.interval=NULL,
    progress=TRUE, logbf = FALSE)
}
\arguments{
  \item{y}{
A vector of observations
}
  \item{iterations}{
Number of MCMC iterations to perform
}
  \item{rscale}{
The prior scale.  A number of preset values can be given as strings; see Details.
}
  \item{null.interval}{
Vector of length two, giving a range of ``negligible'' standardized effect sizes to serve as an interval null hypothesis. Endpoints can be infinite.
}

  \item{progress}{
Report progress with a text progress bar?
}
  \item{t}{
	The t value obtained from a classical t~test (with equal variance assumed, if two-sample)
}
  \item{n1}{
Sample size for first group
}
  \item{n2}{
Sample size for second group. If a one-sample test is desired, \code{n2} should be 0.
}
\item{n}{
Sample size for the single group (Gibbs sampler is for one-sample designs)
}
  \item{prior.cauchy}{
Use a Cauchy prior for the prior on effect size? If false, the prior on effect size is Normal (see Details below) 
}
\item{logbf}{
If \code{TRUE}, return the natural logarithm of the Bayes factor
}
}
\details{
The Bayes factor provided by \code{ttest.Quad} tests the null hypothesis that the mean of a normal population is 0. 
Specifically, the Bayes factor compares two hypotheses: that the standardized effect size is 0, or that the 
standardized effect size is not 0. For one-sample tests, the standardized effect size is mu/sigma; 
for two sample tests, the standardized effect size is (mu_2-mu_1)/sigma.

A noninformative Jeffreys prior is placed on the variance of the normal population, while a Cauchy 
prior is placed on the standardized effect size. The \code{rscale} argument controls the scale of the prior 
distribution, with \code{rscale=1} yielding a standard Cauchy prior. If \code{prior.Cauchy=FALSE}, then a Normal prior is used
instead of a Cauchy prior. See the references below for more details.

For the \code{rscale} argument, several named values are recongized: "medium" corresponds to \eqn{r=\sqrt{2}/2}{r=sqrt(2)/2}; "wide" corresponds to \eqn{r=1}{r=1}.

\code{ttest.Gibbs} samples from the posterior distributions for
one-sample designs. For two (or more) sample designs, see \link{oneWayAOV.Gibbs} or \link{nWayAOV.Gibbs}.
}
\value{
\code{ttest.Quad} returns a scalar giving the Bayes factor in favor of the alternative hypothesis that the effect size is not 0.

\code{ttest.Gibbs} returns a list containing the following:
\item{chains}{A object of type MCMC containing the chains for each parameter}
\item{BF}{An MCMC estimate of the Bayes factor, computed using the Savage-Dickey method (Morey, Rouder, Pratte, and Speckman, 2011)}
\item{BFarea}{If an \code{interval.null} was provided, an MCMC estimate of Bayes factor in favor of the \emph{complement} of the interval null hypothesis.}
}
\references{
Morey, R. D., Rouder, J. N., Pratte, M. S., \& Speckman, P. L. (2011). Using MCMC chain outputs to efficiently estimate Bayes factors. Journal of Mathematical Psychology, 55, 368-378

Morey, R. D. \& Rouder, J. N. (2011). Bayes Factor Approaches for Testing Interval Null Hypotheses. Psychological Methods, 16, 406-419

Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., & Iverson, G. (2009). Bayesian t-tests
for accepting and rejecting the null hypothesis. Psychonomic Bulletin & Review, 16, 752-760

Perception and Cognition Lab (University of Missouri): Bayes factor calculators. \url{http://pcl.missouri.edu/bayesfactor}
}
\author{
Richard D. Morey
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning}{
The default prior scale has changed from 1 to \eqn{\sqrt{2}/2}. The factor of \eqn{\sqrt{2}} is to be consistent with Morey et al. (2011) and Rouder et al. (2012), and the factor of \eqn{1/2} is to better scale the expected effect sizes; the previous scaling put more weight on larger effect sizes. To obtain the same Bayes factors as Rouder et al. (2009), change the prior scale to 1.
}

\seealso{
\code{\link{t.test}}, \code{\link{oneWayAOV.Quad}},\code{\link{oneWayAOV.Gibbs}}, \code{\link{nWayAOV.MC}}, \code{\link{nWayAOV.Gibbs}}
}
\examples{
## Classical example: Student's sleep data
data(sleep)
plot(extra ~ group, data = sleep)

## t.test() gives a t value of -4.0621
t.test(extra ~ group, data = sleep, paired=TRUE)

## Gives a Bayes factor of about 17
## in favor of the alternative hypothesis
ttest.Quad(t = -4.0621, n1 = 10)


## Compute difference scores for one-sample test
diffScores = sleep$extra[11:20] - sleep$extra[1:10]
output = ttest.Gibbs(diffScores)

## Look at the posterior distribution of the mean
plot(output$chains[,1])

## Posterior mean, 95% credible interval
mean(output$chains[,1])
quantile(output$chains[,1],p=c(0.025,0.975))

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{htest}
\keyword{models}