\name{freqCIdelta}
\Rdversion{1.1}
\alias{freqCIdelta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
freqCIdelta(t, N, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
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
function (t, N, alpha = 0.05) 
{
    ow = options("warn")
    options(warn = -1)
    t.hi = optimize(find.ncp.hi, interval = c(-10, 10), tol = .Machine$double.eps^0.25, 
        t = t, alpha = alpha, N = N)$minimum
    t.lo = optimize(find.ncp.lo, interval = c(-10, 10), tol = .Machine$double.eps^0.25, 
        t = t, alpha = alpha, N = N)$minimum
    d.hi = t.hi/sqrt(N)
    d.est = t/sqrt(N)
    d.lo = t.lo/sqrt(N)
    options(ow)
    c(d.lo, d.est, d.hi)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
