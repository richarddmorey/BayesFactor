\name{samp.eta2.t}
\Rdversion{1.1}
\alias{samp.eta2.t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
samp.eta2.t(old, d, t, N, sdmet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{old}{
%%     ~~Describe \code{old} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{sdmet}{
%%     ~~Describe \code{sdmet} here~~
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
function (old, d, t, N, sdmet) 
{
    z = rnorm(1, 0, sdmet)
    cand = old + z
    if (cand <= 0) {
        llr = -Inf
    }
    else {
        llr = proplog.deta2.t(cand, d, t, N) - proplog.deta2.t(old, 
            d, t, N)
    }
    b = rbinom(1, 1, exp(min(llr, 0)))
    b * cand + (1 - b) * old
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
