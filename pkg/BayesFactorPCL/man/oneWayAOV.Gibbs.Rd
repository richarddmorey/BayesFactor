\name{one-way ANOVA}
\alias{oneWayAOV}
\alias{oneWayAOV.Quad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions for Bayesian analysis of one-way designs
}
\description{
\code{oneWayAOV.Quad} computes Bayes factors for balanced one-way designs.
}
\usage{
oneWayAOV.Quad(F, N, J, rscale = "medium", logbf = FALSE)

}
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
	The prior scale. A number of preset values can be given as strings; see Details.
}
\item{logbf}{
If \code{TRUE}, return the natural logarithm of the Bayes factor
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

For the \code{rscale} argument, several named values are recongized: "medium" corresponds to \eqn{r=1/2}{r=1/2}; "wide" corresponds to \eqn{r=\sqrt{2}/2}{r=sqrt(2)/2}.


}
\value{
\code{oneWayAOV.Quad} returns a scalar Bayes factor in favor of the alternative hypothesis.
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
\code{oneWayAOV.Quad} should only be used with F values obtained from balanced designs. For unbalanced designs, use \code{\link{nWayAOV.MC}}.
}
\seealso{
\code{\link{aov}}, \code{\link{ttest.Quad}}, \code{\link{nWayAOV.MC}}
}

\examples{
## Example data "InsectSprays" - see ?InsectSprays
require(stats); require(graphics)
boxplot(count ~ spray, data = InsectSprays,
        xlab = "Type of spray", ylab = "Insect count",
        main = "InsectSprays data", varwidth = TRUE, col = "lightgray")

## Classical analysis (with transformation)
classical <- aov(sqrt(count) ~ spray, data = InsectSprays)
plot(classical)
summary(classical)

## Bayes factor (a very large number)
Fvalue <- anova(classical)$"F value"[1]
oneWayAOV.Quad(Fvalue, N=12, J=6)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
