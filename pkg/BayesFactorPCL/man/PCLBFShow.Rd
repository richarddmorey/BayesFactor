\name{PCLBFShow}
\Rdversion{1.1}
\alias{PCLBFShow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PCLBFShow(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object) 
{
    priorType = ifelse(object@prior.cauchy, "Cauchy", "Normal")
    priorType = paste(priorType, " (scale of ", object@prior.scale, 
        ") ", sep = "")
    intType = ifelse(object@MCMC, paste("MCMC (", object@MCMC.iters, 
        " iterations)", sep = ""), "Adaptive Quadrature")
    mcmc95 = ifelse(object@MCMC, paste("*The MCMC 95\% CI on the NOH Bayes Factor is (", 
        zapsmall(object@MCMC.95CI[1]), ",", zapsmall(object@MCMC.95CI[2]), 
        ")\n", sep = ""), "")
    display = paste("Bayesian Analysis\n-----------------\n", 
        "Prior type         : ", priorType, "\n", "JZS Bayes Factor   : ", 
        object@JZS, "\n", "NOH Bayes Factor   : ", object@NOH, 
        " (null region is ", -object@null.d, ",", object@null.d, 
        ")\n", mcmc95, "Hybrid Bayes Factor: ", object@hybrid, 
        " (pi0 is ", object@pi0, "; null region is ", -object@null.d, 
        ",", object@null.d, ")\n", "\n\nIntegration Type: ", 
        intType, "\n\n", "Frequentist Analysis\n--------------------\n", 
        100 * object@freq.conf, "\% CI on d: (", zapsmall(object@freq.CI[1]), 
        ",", zapsmall(object@freq.CI[3]), ")\n", "p value    : ", 
        zapsmall(object@freq.p), "\n\n", sep = "")
    cat(noquote(display))
    return(NULL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
