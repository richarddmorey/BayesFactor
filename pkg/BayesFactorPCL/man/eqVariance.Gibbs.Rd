\name{eqVariance.Gibbs}
\alias{eqVariance.Gibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Bayes factor test of whether variances of multiple groups are equal 
}
\description{
This function computes a Bayes factor hypothesis test for whether the variances of multiple groups are equal, 
and also samples from the posterior distribution of the 
}
\usage{
eqVariance.Gibbs(y, iterations = 1000, lambda = 0.1, M2scale = 1, sig2.metrop.sd = 1, tau.metrop.sd = 1, M2.metrop.scale = 1, newtonSteps = 6, progress = TRUE, whichModel = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	A matrix of data, with groups as columns, and each row is an observation. 
	The total number of rows should be the maximum number of observations in a group. If some groups have fewer 
	observations than others, then those columns should be padded with NAs.
}
  \item{iterations}{
Number of MCMC iterations to perform
}
  \item{lambda}{
	Only used if \code{whichModel=1}. The rate of the exponential prior distribution on $\tau$.
  }
  \item{M2scale}{
	Only used if \code{whichModel=2}. Prior scale for $g$ under Model~2.
}
  \item{sig2.metrop.sd}{
	Only used if \code{whichModel=1}. Scale for the Metropolis-Hastings sampling of $\sigma^2$.
}
  \item{tau.metrop.sd}{
Only used if \code{whichModel=1}. Scale for the Metropolis-Hastings sampling of $\tau$.
}
  \item{M2.metrop.scale}{
	Only used if \code{whichModel=2}. Scale for the Metropolis-Hastings sampling of $g$.
}
  \item{newtonSteps}{
	Only used if \code{whichModel=2}. The number of Newton-Raphson steps to use in finding the mode of the full conditional posterior distribution of $g_j$. See Details below.
}
  \item{progress}{
Report progress with a text progress bar?
 }
  \item{whichModel}{
	Which prior structure on the group variance scales to use? \code{whichModel=1} uses a conjugate gamma 
	distribution, but is not recommended. \code{whichModel=2} uses a log-normal distribution, and is recommended. See Details below.
}
}
\details{
This function computes a Bayes factor test for equal variances for an arbitrary number of independent groups. All observations are 
assumed normally distributed; each group is allowed to have its own mean $\mu_j$ and variance $\sigma^2_j = g_j\sigma^2$.
Noninformative Jeffreys priors are placed on $\mu_j$ and $\sigma^2$.

Under Model~1, the prior distribution on $g$ is a gamma with rate and shape both equal to $\tau$. 
An exponential prior with rate \code{lambda} is placed on $\tau$. This model is not recommended because the gamma distribution lacks multiplicative 
symmetry.

Under Model~2, the prior distribution of the logarithm of $g$ is multivariate t(1) distribution, with scale \code{M2scale}. In estimating the Bayes factor,
the IWMDE method of Chen (1994) was implemented. The argument \code{newtonSteps} controls the quality of the approximation used in the IWMDE method:
higher numbers of steps lead to more efficient estimates of the Bayes factor. See references for more details.

}
\value{
\item{chains}{
An object of class mcmc, containing MCMC chains for all parameters
}

\item{bayesFactor}{
The Bayes factor in favor of the null hypothesis that all group variances are equal
}

\item{acc.rates}{
Acceptance rates for the Metropolis-Hastings steps
}
}
\references{
\item{}{Morey, R. D., Rouder, J. N., Pratte, M. S., \& Speckman, P. L. (submitted). Conditional Marginal Density 
Estimation to Improve Bayes Factor Computation.}
}
\author{
Richard D. Morey
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LeveneTest}} in the \code{car} package
}
\examples{
## Use the singer heights data to 
## test the hypothesis that males and female
## singers have different population variances
data(singer)

## Format for the equal variance test
y = matrix(NA,128,2)
## Females
y[,1] = singer[1:128,1]
## Males
y[1:107,2] = singer[-(1:128),1]

## Make boxplot
boxplot(y,names=c("Females","Males"))

## Analysis
output = eqVariance.Gibbs(y,whichModel=2,M2.metrop.scale=2)

## Equal variance Bayes factor
## Should be about 1.8 in favor of null
output[[2]]

## Levene's test for comparison
require(car) # For Levene's test
singer$sex = as.factor(c(rep('F',128),rep('M',107)))

## p value is 0.1855
leveneTest(height~sex,data=singer)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
