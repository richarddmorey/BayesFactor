\name{ols.value.bf}
\Rdversion{1.1}
\alias{ols.value.bf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ols.value.bf(F, n, v = 1, prior.cauchy = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{F}{
%%     ~~Describe \code{F} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{v}{
%%     ~~Describe \code{v} here~~
}
  \item{prior.cauchy}{
%%     ~~Describe \code{prior.cauchy} here~~
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
function (F, n, v = 1, prior.cauchy = T) 
{
    Rsqr = F/(F + n - 2)
    bfinv = ifelse(prior.cauchy, integrate(ols.g, lower = 0, 
        upper = Inf, Rsqr = Rsqr, n = n, v = v)$value, exp(((n - 
        2)/2) * log(1 + n * v) + (-(n - 1)/2) * log(1 + (1 - 
        Rsqr) * n * v)))
    return(1/bfinv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
