\name{find.crit.t.areabf}
\Rdversion{1.1}
\alias{find.crit.t.areabf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find.crit.t.areabf(bf = 1/10, N, null.d = 0.2, df.prior = 1, bounds = c(0, 20), r.scale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bf}{
%%     ~~Describe \code{bf} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{null.d}{
%%     ~~Describe \code{null.d} here~~
}
  \item{df.prior}{
%%     ~~Describe \code{df.prior} here~~
}
  \item{bounds}{
%%     ~~Describe \code{bounds} here~~
}
  \item{r.scale}{
%%     ~~Describe \code{r.scale} here~~
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
function (bf = 1/10, N, null.d = 0.2, df.prior = 1, bounds = c(0, 
    20), r.scale = 1) 
optimize(critt.work, interval = bounds, N = N, null.d = null.d, 
    df.prior = df.prior, r.scale = r.scale)$minimum
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
