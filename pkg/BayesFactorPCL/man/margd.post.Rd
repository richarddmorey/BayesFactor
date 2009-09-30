\name{margd.post}
\Rdversion{1.1}
\alias{margd.post}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
margd.post(delta, t, N, df.prior = 1, eta2.lo = 0, eta2.up = Inf, r.scale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{df.prior}{
%%     ~~Describe \code{df.prior} here~~
}
  \item{eta2.lo}{
%%     ~~Describe \code{eta2.lo} here~~
}
  \item{eta2.up}{
%%     ~~Describe \code{eta2.up} here~~
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
function (delta, t, N, df.prior = 1, eta2.lo = 0, eta2.up = Inf, 
    r.scale = 1) 
integrate(vect.jnt.post, eta2.lo, eta2.up, delta = delta, t = t, 
    N = N, df.prior = df.prior, r.scale = r.scale, rel.tol = max(50 * 
        .Machine$double.eps, 5e-29))$value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
