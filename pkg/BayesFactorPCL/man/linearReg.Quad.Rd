\name{linearReg.Quad}
\alias{linearRed.Quad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to compute the Bayes factors regression designs
}
\description{
This function computes the Bayes factor for multiple linear regression, testing whether slopes are zero or not.
}
\usage{
linearReg.Quad(N,p,R2)
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

}
\details{
This Bayes factor tests the hypothesis that all slopes are 0. 

}
\value{
A scalar giving the Bayes factor in favor of the alternative hypothesis at least one slope is nonzero.
}
\references{
\item{}{Liang, F. and Paulo, R. and Molina, G. and Clyde, M. A. and Berger, J. O. (2008). Mixtures of g-priors for Bayesian Variable Selection.
Journal of the American Statistical Association, 103, pp. 410-423}

\item{}{Rouder, J. N.  and Morey, R. D. (submitted). Bayesian testing in regression.}

\item{}{Perception and Cognition Lab (University of Missouri): Bayes factor calculators. \link{http://pcl.missouri.edu/bayesfactor}}
}
\author{
Jeffrey N. Rouder
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lm}}
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
## the data strongly favor the alternative
linearReg.Quad(30,1,0.6813)


}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
