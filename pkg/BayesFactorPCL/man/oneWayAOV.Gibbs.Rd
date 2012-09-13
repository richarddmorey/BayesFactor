\name{oneWayAOV.Gibbs}
\alias{oneWayAOV.Gibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to sample from the posterior distribution of parameters in one-way designs 
}
\description{
This function samples from the posterior distributions of parameters for 
one-way designs, using Gibbs sampling.  See details about priors in the details.
}
\usage{
oneWayAOV.Gibbs(y, iterations = 10000, rscale = 1, progress = TRUE, gibi=NULL)
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
}
\details{
For details about the priors, see help for the function \link{oneWayAOV.Quad}.
}
\value{
\item{chains}{A object of type MCMC containing the chains for each parameter}
\item{BF}{An MCMC estimate of the Bayes factor, computed using the Savage-Dickey method (Morey, Rouder, Pratte, and Speckman, submitted)}
}

\references{
Morey, R. D., Rouder, J. N., Pratte, M. S., \& Speckman, P. L. (submitted). Conditional Marginal Density 
Estimation to Improve Bayes Factor Computation.

Perception and Cognition Lab (University of Missouri): Bayes factor calculators. \url{http://pcl.missouri.edu/bayesfactor}
}
\author{
Richard D. Morey
}
\note{
For more accurate estimation of the Bayes factor by Gaussian quadrature, see \code{\link{oneWayAOV.Quad}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{aov}}, \code{\link{ttest.Quad}}, \code{\link{oneWayAOV.Quad}}
}

\examples{
## Classical example (from t.test help):
## Student's sleep data
data(sleep)
plot(extra ~ group, data = sleep)

## Format for oneWayAOV.Gibbs
y = matrix(sleep[,1],10,2) 
output = oneWayAOV.Gibbs(y)

## Bayes factor should be about 1.02
output[[2]]

## Plot posteriors of grand mean and effects
plot(output[[1]][,1:3])

## Standard ANOVA analysis
## yields p value of about 0.08
summary(aov(extra ~ group, data = sleep))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
