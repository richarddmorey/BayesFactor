\name{allNways}
\alias{allNways}
\alias{nWayAOV.MC}
\alias{nWayAOV.Gibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to compute the Bayes factors for general ANOVA designs
}
\description{
These functions computes Bayes factors which enable testing main effects and interactions in general ANOVA designs. \code{allNways} will compute Bayes factors for testing all fixed effects in mixed ANOVA models against the null hypothesis that all fixed effects are 0, and \code{nWayAOV.MC} computes the Bayes factor for an arbitrary design matrix.
}
\usage{
allNways(y, dataFixed = NULL, dataRandom = NULL, iterations = 1000, 
    samples = FALSE, only.top = FALSE, ...) 
    
nWayAOV.MC(y, X, struc, iterations = 10000, rscale = 1, progress = FALSE, 
    samples = FALSE, gsamples = FALSE, gibi = NULL)
    
nWayAOV.Gibbs(y, X, struc, iterations = 10000, rscale = 1, progress = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	A vector of observations	
}
  \item{dataFixed}{
A data frame containing the fixed factors, one column per factor. Columns must be factors
}
  \item{dataRandom}{
A data frame containing the fixed factors, one column per factor. Columns must be factors
}
  \item{iterations}{
Number of Monte Carlo samples to perform (see details)
}
  \item{samples}{
If \code{TRUE}, return the Monte Carlo samples used in estimating the integral (default is \code{FALSE})
}
  \item{only.top}{
If \code{TRUE}, only the most complex models will be tested, cutting down the computation time (see details)
}
\item{X}{
A design matrix whose number of rows match \code{length(y)}
}
\item{struc}{
Vector grouping the columns of \code{X} (see details)
}
\item{rscale}{
Either a scalar or a vector of the same length as \code{struc}, of prior scales (see details)
}
\item{progress}{
If \code{TRUE}, text progress bar is displayed
}
\item{gsamples}{
If \code{TRUE}, return all the Monte Carlo samples of the \code{g} parameter drawn (default is \code{FALSE}; see details)
}
\item{gibi}{
Interface for a future graphical user interface (not intended for use by end users)
}
\item{...}{
Additional arguments to pass to \code{nWayAOV.MC}
}
}
\details{

}
\value{
If \code{samples} and \code{gsamples} are both \code{FALSE}, then \code{allNways} and \code{nWayAOV.MC} return a vector of log Bayes factor or a single log10 Bayes factors, respectively. 

If either \code{samples} or \code{gsamples} are \code{TRUE}, the \code{nWayAOV.MC} will return lists whose first element is as above, and the subsequent elements contains a vector with the samples. \code{allNways} will return a list whose first element is a vector as above, and whose second element is a list containing the samples from the respective models.
}

\references{
Rouder, J. N., Morey, R. D., Speckman, P. L., Province, J. M., (in press) Default Bayes Factors for ANOVA Designs. Journal of Mathematical Psychology

Perception and Cognition Lab (University of Missouri): Bayes factor calculators. \url{http://pcl.missouri.edu/bayesfactor}
}
\author{
Richard D. Morey
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{
The function \code{allNways} will compute Bayes factors for all possible combinations of fixed factors and interactions, against the null hypothesis that \emph{all} effects are 0. The total number of tests computed will be \eqn{2^{2^K - 1}}{2^(2^K - 1)} for \eqn{K} fixed factors. This number increases very quickly with the number of factors. For instance, for a five-way ANOVA, the total number of tests exceeds two billion. Even though each test takes a fraction of a second, the time taken for all tests could exceed your lifetime.

It is possible to reduce the number of models tested by only testing the most complex model and every restriction that can be formed by removing one factor or interaction using the \code{only.top} argument. This reduces the number of tests to \eqn{2^K-1}, which is more manageable. The Bayes factor for each restriction against the most complex model can be interpreted as a test of the removed factor/interaction. If other tests are required besides those provided by \code{only.top}, they can be obtained via the function \code{nWayAOV.MC}.

}
\seealso{
\code{\link{aov}}, \code{\link{ttest.Quad}}, \code{\link{oneWayAOV.Quad}},\code{\link{oneWayAOV.Gibbs}}
}
\examples{
## Classical example, taken from t.test() example
## Student's sleep data
data(sleep)
plot(extra ~ group, data = sleep)

## traditional ANOVA gives a p value of 0.00283
summary(aov(extra ~ group + Error(ID), data = sleep))

logbf <- allNways(y = sleep$extra,
          dataFixed = data.frame(sleep$group), 	   
          dataRandom = data.frame(sleep$ID), 
          iterations=10000)

## Gives a Bayes factor of about 14
## in favor of the alternative hypothesis
exp(logbf)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
