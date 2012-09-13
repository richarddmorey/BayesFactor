\name{oneWayAOV.Gibbs}
\alias{oneWayAOV.Gibbs}
\alias{oneWayAOV.Quad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions for Bayesian analysis of one-way designs
}
\description{
\code{oneWayAOV.Gibbs} samples from the posterior distributions of parameters for 
one-way designs, using Gibbs sampling; \code{oneWayAOV.Quad} computes Bayes factors for balanced designs.
}
\usage{
oneWayAOV.Gibbs(y, iterations = 10000, rscale = 1, progress = TRUE, gibi=NULL)

oneWayAOV.Quad(F, N, J, rscale = 1)
}
\arguments{
  \item{y}{
	A matrix of data, with groups as columns, and each row is an observation. 
	The total number of rows should be the maximum number of observations in a group. If some groups have fewer 
	observations than others, then those columns should be padded with NAs.
}
  \item{iterations}{
Number of MCMC iterations to perform
}
  \item{rscale}{
	The prior scale (see Details below)
}
  \item{progress}{
Report progress with a text progress bar?
}
\item{gibi}{
Interface for a future graphical user interface (not intended for use by end users)
}
\item{F}{
	The F value from the traditional ANOVA analysis
}
  \item{N}{
	The number of observations per group
}
  \item{J}{
	The number of groups
}
}
\details{
This Bayes factor tests the null hypothesis that the mean of a multivariate normal population is 0. 
Specifically, the Bayes factor compares two hypotheses: that the standardized effect size is 0, or that the 
standardized effect size is not 0. Here, the standardized effect size is mu/sigma where mu is the vector of true 
group means and sigma is the error variance. 

A noninformative Jeffreys prior is placed on the error variance, while a multivariate t(1) prior is placed on
standardized effect size (\eqn{\mu_j/\sigma}{mu_j/sigma}, for the \eqn{j}{j}th group) . The \code{rscale} argument controls the scale of the prior 
distribution, with \code{rscale=1} yielding a standard multivariate t(1). See the references below for more details.

}
\value{
\code{oneWayAOV.Gibbs} returns a list containing two elements:
\item{chains}{A object of type MCMC containing the chains for each parameter}
\item{BF}{An MCMC estimate of the Bayes factor, computed using the Savage-Dickey method (Morey, Rouder, Pratte, and Speckman, 2011)}

\code{oneWayAOV.Quad} returns a scalar Bayes factor in favor of the null hypothesis.
}

\references{
Morey, R. D., Rouder, J. N., Pratte, M. S., \& Speckman, P. L. (2011). Using MCMC chain outputs to efficiently estimate Bayes factors. Journal of Mathematical Psychology, 55, 368-378

Perception and Cognition Lab (University of Missouri): Bayes factor calculators. \url{http://pcl.missouri.edu/bayesfactor}
}
\author{
Richard D. Morey
}


%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{
\code{oneWayAOV.Quad} should only be used with F values obtained from balanced designs. For unbalanced designs, use \code{oneWayAOV.Gibbs} or \code{\link{nWayAOV.MC}}.
}
\seealso{
\code{\link{aov}}, \code{\link{ttest.Quad}}, \code{\link{nWayAOV.MC}}
}

\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
