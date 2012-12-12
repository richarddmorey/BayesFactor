\name{Linear Regression}
\alias{linearReg.Quad}
\alias{linearReg.Gibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A functions to compute the Bayes factors and sample from posteriors in regression designs
}
\description{
These functions compute the Bayes factor for multiple linear regression, testing whether slopes are zero or not, or sample from the corresponding posteriors.
}
\usage{
linearReg.Quad(N,p,R2, rscale=1, logbf=FALSE)

linearReg.Gibbs(y, covariates, iterations = 10000, rscale = 1, 
               progress = TRUE, gibi=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
	Total number of observations
}
  \item{p}{
	Total number of predictors
}
\item{R2}{
  The coefficient of detemination for the full regression model containing all \code{p} predictors
}
\item{rscale}{
  Scale of the Cauchy prior on standardized effect size  
}
\item{logbf}{
If \code{TRUE}, return the natural logarithm of the Bayes factor
}
\item{y}{
  A vector of observations	
}
\item{covariates}{
A data frame containing the covariates, one column per covariate. 
}
\item{iterations}{
Number of Gibbs sampler iterations to perform
}
\item{progress}{
If \code{TRUE}, text progress bar is displayed
}
\item{gibi}{
Interface for a future graphical user interface (not intended for use by end users)
}
}
\details{
\code{linearReg.Quad} computes Bayes factors that test the hypothesis that all slopes are 0 against the alternative that all slopes are nonzero. \code{linearReg.Quad} uses the Gibbs sampler to sample from the posterior distribution assuming that all slopes are non-zero. 

The vector of observations \eqn{y} is assumed to be distributed as
\deqn{y ~ Normal(\alpha 1 + X\beta, \sigma^2 I).}
The joint prior on \eqn{\alpha,\sigma^2} is proportional to \eqn{1/\sigma^2}, the prior on \eqn{\beta} is
\deqn{\beta ~ Normal(0, N g \sigma^2(X'X)^{-1}).}
where \eqn{g ~ InverseGamma(1/2,r/2)}. See Liang et al. (2008) section 3 for details.
}
\value{
\code{linearReg.Quad} returns a scalar giving the Bayes factor in favor of the alternative hypothesis at least one slope is nonzero. 

\code{linearReg.Gibbs} returns an \code{mcmc} object containing the samples from the Gibbs sampler. 
}
\references{
Liang, F. and Paulo, R. and Molina, G. and Clyde, M. A. and Berger, J. O. (2008). Mixtures of g-priors for Bayesian Variable Selection.
Journal of the American Statistical Association, 103, pp. 410-423

Rouder, J. N.  and Morey, R. D. (submitted). Bayesian testing in regression.

Perception and Cognition Lab (University of Missouri): Bayes factor calculators. \url{http://pcl.missouri.edu/bayesfactor}
}
\author{
Jeffrey N. Rouder, Richard D. Morey
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lm}}, \code{\link{mcmc}}
}
\examples{
## Use attitude data set
data(attitude)

## Scatterplot
lm1 = lm(rating~complaints,data=attitude)
plot(attitude$complaints,attitude$rating)
abline(lm1)

## Traditional analysis
## p value is highly significant
summary(lm1)

## Bayes factor
## The Bayes factor is 792240.6;
## the data strongly favor hypothesis that 
## the slope is not 0.
linearReg.Quad(30,1,0.6813)


}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
