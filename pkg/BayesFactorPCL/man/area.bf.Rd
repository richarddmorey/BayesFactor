\name{area.bf}
\Rdversion{1.1}
\alias{area.bf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
area.bf(t, N1, N2 = NULL, df.prior = 1, null.d = 0.2, MCMC = FALSE, r.scale = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{N1}{
%%     ~~Describe \code{N1} here~~
}
  \item{N2}{
%%     ~~Describe \code{N2} here~~
}
  \item{df.prior}{
%%     ~~Describe \code{df.prior} here~~
}
  \item{null.d}{
%%     ~~Describe \code{null.d} here~~
}
  \item{MCMC}{
%%     ~~Describe \code{MCMC} here~~
}
  \item{r.scale}{
%%     ~~Describe \code{r.scale} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (t, N1, N2 = NULL, df.prior = 1, null.d = 0.2, MCMC = FALSE, 
    r.scale = 1, ...) 
{
    if (MCMC) {
        out = area.bf.mcmc(t, N1, df.prior, null.d, , r.scale = r.scale, 
            ...)
    }
    else {
        out = area.bf.integrate(t, N1, df.prior, null.d, r.scale = r.scale, 
            ...)
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
