\name{nWayAOV.MC}
\alias{nWayAOV.MC}
\alias{nWayAOV}
\alias{nWayAOV.Gibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to compute Bayes factors for a specified ANOVA model defined by a design matrix.
}
\description{
These functions computes a single Bayes factor, or estimate parameters, for a specified ANOVA model to a null ANOVA model defined by the addition of a grand mean and error terms.}
\usage{
    
nWayAOV.MC(y, X, struc, iterations = 10000, rscale = 1, progress = FALSE, 
    samples = FALSE, gsamples = FALSE, gibi = NULL, logbf = FALSE)
    
nWayAOV.Gibbs(y, X=NULL, struc=NULL, dataFixed=NULL, dataRandom=NULL, iterations = 10000, rscale = 1, progress = TRUE, unreduce=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{
	A vector of observations	
}
\item{X}{
A design matrix whose number of rows match \code{length(y)}.
}

\item{struc}{
Vector grouping the columns of \code{X} (see details).
}
  \item{dataFixed}{
An (optional) data frame containing the fixed factors, one column per factor. Columns will be converted to factors if they are not factors already. If \code{dataFixed} or \code{dataRandom} are specified, \code{X} and \code{struc} for the corresponding full model will be constructed
}
  \item{dataRandom}{
An (optional) data frame containing the random factors, one column per factor. Columns will be converted to factors if they are not factors already.
}

\item{iterations}{
Number of Monte Carlo samples used to estimate Bayes factor.
}

\item{rscale}{
A scalar or a vector (same length as \code{struc}) of prior scale(s) (see details).
}

\item{progress}{
If \code{TRUE}, text progress bar is displayed
}
\item{unreduce}{
  \code{TRUE}, project the \code{K-1} effects back into \code{K} effects that sum to 0 (if \code{X} is constructed using \code{dataFixed})
}
  \item{samples}{
If \code{TRUE}, return the Monte Carlo samples used in estimating the integral (default is \code{FALSE})
}
 

\item{gsamples}{
If \code{TRUE}, return all the Monte Carlo samples of the \code{g} parameter drawn (default is \code{FALSE}; see details)
}
\item{gibi}{
Interface for a future graphical user interface (not intended for use by end users)
}
\item{logbf}{
If \code{TRUE}, return the natural logarithm of the Bayes factor
}

}
\details{
The models, priors, and methods of computation are provided in Rouder et al. (2012).  	
The linear model for a vector of \eqn{N} observations, denoted by \eqn{y}, is 
\deqn{ y = \mu + X \theta +\epsilon,}
where \eqn{\theta} is a vector of \eqn{p} effects, \eqn{X} is an \eqn{N}-by-\eqn{p} design matrix, and \eqn{\epsilon} is a vector of \eqn{N} zero-centered noise terms with variance \eqn{\sigma^2}.    Zellner and Siow (1980) inspired g-priors are placed on effects, but with a possibly separate g-prior parameter for each effect:
\deqn{\theta~N(0,\sigma^2G),} 
where \eqn{G} is a \eqn{p}-by-\eqn{p} diagonal matrix.  

Following Gelman (2005), effects may grouped into blocks where effects within a single block are exchangable and those across blocks are not.  These blocks correspond to the main effect of a factor or the interaction between factors.  Blocking is implemented by placing constraints on diag(G) as follows.  Consider a design with two factors each with two levels.  There are eight effects (2 main effects for the two levels and 4 interaction effects) which should be grouped into three blocks (one for each main effect, and one for the interaction).  Exchangability within blocks but not across it is given by \eqn{diag(G)=(g_1, g_1, g_2, g_2, g_3, g_3, g_3, g_3, g_3)} where the first two elements are for one main effect, the next two are for the other, and the last four are for the interaction (assuming this order is coded in the design matrix \eqn{X}).  \code{struct} codes how the diagonal of \eqn{G} should be grouped, with each entry indicating how many levels are in a block. Note \code{sum(struct)} must equal the number of non-constant columns of \code{X}.  

A Jeffreys prior is placed on \eqn{\mu} and \eqn{\sigma^2}.  Independent scaled inverse-chi-square priors with one degree of freedom are placed on \eqn{g_1,\ldots,g_r}, where \eqn{r} is the number of blocks.  The square-root of the scale for \eqn{g}'s is given in \code{rscale}.

The null model is given by \deqn{y_i = \mu + \epsilon_i} for \eqn{i=1,...,N}.

\code{nWayAOV.Gibbs} returns posterior samples of the parameters, which are computed using Markov-chain Monte Carlo sampling.  
}

\value{
With default settings, the Bayes factor for the model relative to the null is returned.  If either \code{samples} or \code{gsamples} are \code{TRUE}, then a list whose first element is as above, and the subsequent elements contain a vector containing the samples.
}

\references{
Gelman, A. (2005) Analysis of Variance---why it is more important than ever.  Annals of Statistics, 33, pp. 1-53. 	
	
Rouder, J. N., Morey, R. D., Speckman, P. L., Province, J. M., (2012) Default Bayes Factors for ANOVA Designs. Journal of Mathematical Psychology.  56.  p. 356-374.

Zellner, A. and Siow, A., (1980) Posterior Odds Ratios for Selected Regression Hypotheses.  In Bayesian Statistics: Proceedings of the First Interanational Meeting held in Valencia (Spain).  Bernardo, J. M., Lindley, D. V., and Smith A. F. M. (eds), pp. 585-603.  University of Valencia.
}
\author{
Richard D. Morey, Jeffrey N. Rouder
}
\note{\code{allNways} is a front-end function for computing multiple Bayes factors in ANOVA designs.  For each submodel, \code{allNways} constructs objects \code{X}, \code{struct}, and \code{rscale} and calls \code{nWayAOV.MC}.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{allNways}}, 	
\code{\link{aov}}, \code{\link{ttest.Quad}}, \code{\link{oneWayAOV.Quad}},\code{\link{oneWayAOV.Gibbs}}
}
\examples{
## Classical example, taken from t.test() example
## Student's sleep data
data(sleep)
plot(extra ~ group, data = sleep)

## traditional ANOVA gives a p value of 0.00283
summary(aov(extra ~ group + Error(ID/group), data = sleep))

## Build design matrix
group.column <- rep(1/c(-sqrt(2),sqrt(2)),each=10)
subject.matrix <- model.matrix(~sleep$ID - 1,data=sleep$ID)
X <- cbind(1, group.column, subject.matrix)

## Bayes factor of full model against grand-mean only model
bf.full <- nWayAOV.MC(y = sleep$extra, X = X, struc = c(1,10), rscale=c(.5,1))

## Sample from posterior distribution of full model
## And plot estimated effect of drug
samples.full <- nWayAOV.Gibbs(y = sleep$extra, X = X, struc = c(1,10), rscale = c(0.5,1), progress = FALSE)
## transform effect so that it is interpretable as difference in means, and plot
plot(samples.full[,2] * sqrt(2))

\dontrun{
## Alternative (easier) way to sample from the 
## posterior distribution of the full model
samples.full2 = nWayAOV.Gibbs(y = sleep$extra,
        dataFixed = data.frame(sleep$group),      
        dataRandom = data.frame(sleep$ID),
        rscale = c(0.5,1), progress = FALSE)
}

## Remove subject effect
X <- cbind(1, subject.matrix)
## Bayes factor of subject-only model against grand-mean only model
bf.reduced <- nWayAOV.MC(y = sleep$extra, X = X, struc = 10, rscale=1)

## Bayes factor test of drug effect
bf.full / bf.reduced

## Compare with allNways result (should be about the same)
allNways(y = sleep$extra,
      dataFixed = data.frame(sleep$group), 	   
      dataRandom = data.frame(sleep$ID), progress=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
