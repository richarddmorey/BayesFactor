\name{oneWayAOV.Quad}
\alias{oneWayAOV.Quad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to compute the Bayes factors for one-way designs
}
\description{
This function computes the Bayes factor corresponding to a balanced one-way ANOVA, using Gaussian quadrature integration. 
See details about priors in the details.
}
\usage{
oneWayAOV.Quad(F, N, J, rscale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{F}{
	The F value from the traditional ANOVA analysis
  }
  \item{N}{
	The number of observations per group
}
  \item{J}{
	The number of groups
}
  \item{rscale}{
	The prior scale (see Details below)
}
}
\details{
This Bayes factor tests the null hypothesis that the mean of a multivariate normal population is 0. 
Specifically, the Bayes factor compares two hypotheses: that the standardized effect size is 0, or that the 
standardized effect size is not 0. Here, the standardized effect size is $\mu/\sigma$ where $\mu$ is the vector of true 
group means and $\sigma$ is the error variance. 

A noninformative Jeffreys prior is placed on the error variance, while a multivariate t(1) prior is placed on
standardized effect size. The \code{rscale} argument controls the scale of the prior 
distribution, with \code{rscale=1} yielding a standard multivariate t(1). See the references below for more details.
}
\value{
A scalar giving the Bayes factor in favor of the null hypothesis that the standardized effect size is 0.
}
\references{
\item{}{Morey, R. D., Rouder, J. N., Pratte, M. S., \& Speckman, P. L. (submitted). Conditional Marginal Density 
Estimation to Improve Bayes Factor Computation.}

\item{}{Perception and Cognition Lab (University of Missouri): Bayes factor calculators. \link{http://pcl.missouri.edu/bayesfactor}}
}
\author{
Richard D. Morey
}
\note{
For unbalanced designs, and samples from the posterior distributions of the parameters, see \code{\link{oneWayAOV.Gibbs}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{aov}}, \code{\link{ttest.Gibbs}}, \code{\link{oneWayAOV.Gibbs}}
}
\examples{
## Classical example (from t.test help):
## Student's sleep data
data(sleep)
plot(extra ~ group, data = sleep)

## Standard ANOVA analysis
## yields an F of about 3.46 
## and p value of about 0.08
summary(aov(extra ~ group, data = sleep))

## Yields a Bayes factor of about 1.02
oneWayAOV.Quad(3.46,10,2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
