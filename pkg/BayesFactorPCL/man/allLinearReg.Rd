\name{allLinearReg}
\alias{allLinearReg}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the Bayes factors in a regression design
}
\description{
This function computes Bayes factors for all models in a regression design
}
\usage{
allLinearReg(y, covariates, extraInfo = FALSE, progress = FALSE, 
             rscale = 1, logbf = FALSE, ...)
    
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	A vector of observations	
}
\item{covariates}{
A data frame containing the covariates, one column per covariate. 
}
\item{extraInfo}{
  If \code{TRUE}, return a data frame containing useful information. If \code{FALSE}, only the Bayes factors are returned. 
}
\item{progress}{
If \code{TRUE}, a text progress bar is displayed. A callback function may also be supplied.
}
\item{rscale}{The prior scale on standardized effect size.
}
\item{logbf}{
If \code{TRUE}, return the natural logarithm of the Bayes factor
}
\item{...}{
Additional arguments to pass to \code{\link{linearReg.Quad}} (currently unused)
}
}

\details{
The models, priors, and methods of computation are provided in Liang et al. (2008). See \code{\link{linearReg.Quad}} for more details. 
}

\value{
If \code{extraInfo} is \code{FALSE}, a vector of Bayes factors, each measuring the support for a model relative to the null model. 

If \code{extraInfo} is \code{TRUE}, a data frame whose three columns are the Bayes factors, R-squared, and number of covariates respectively. Each row is a different model, as indicated by the corresponding row name.
}

\references{
Liang, F. and Paulo, R. and Molina, G. and Clyde, M. A. and Berger, J. O. (2008). Mixtures of g-priors for Bayesian Variable Selection.
Journal of the American Statistical Association, 103, pp. 410-423

Rouder, J. N.  and Morey, R. D. (in press). Bayesian testing in regression. Multivariate Behavioral Research.

Zellner, A. and Siow, A., (1980) Posterior Odds Ratios for Selected Regression Hypotheses.  In Bayesian Statistics: Proceedings of the First Interanational Meeting held in Valencia (Spain).  Bernardo, J. M., Lindley, D. V., and Smith A. F. M. (eds), pp. 585-603.  University of Valencia.

}
\author{
Richard D. Morey, Jeffrey N. Rouder
}



%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{
The function \code{allLinearReg} will compute Bayes factors for all combinations of covariates. The total number of models tested will be \eqn{2^K - 1}{2^K - 1} for \eqn{K} covariates. For example, with 16 covariates, over 65,000 models will be tested. Even though each test takes a fraction of a second, for large numbers of covariates the function could take a long time time to finish. Use the \code{progress} argument to keep an eye on the function's progress.
}

\seealso{
\code{\link{linearReg.Quad}}
}
\examples{
## See help(attitude) for details about the data set
data(attitude)

## Classical regression
summary(fm1 <- lm(rating ~ ., data = attitude))

## Compute Bayes factors for all regression models
allLinearReg(y = attitude$rating, covariates = attitude[,-1])
## Best model is 'complaints' only

## Output extra info, and compute all Bayes factors against full model
output <- allLinearReg(y = attitude$rating, covariates = attitude[,-1],
    extraInfo = TRUE)
output$bf <- output$bf / output$bf[which.max(output$p)]

## What are the top models?
head(output)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
