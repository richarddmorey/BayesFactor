\name{hybrid.bf}
\Rdversion{1.1}
\alias{hybrid.bf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hybrid.bf(t, N, pi0 = 0.5, prior.cauchy = TRUE, null.d = 0.2, MCMC = FALSE, scale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{pi0}{
%%     ~~Describe \code{pi0} here~~
}
  \item{prior.cauchy}{
%%     ~~Describe \code{prior.cauchy} here~~
}
  \item{null.d}{
%%     ~~Describe \code{null.d} here~~
}
  \item{MCMC}{
%%     ~~Describe \code{MCMC} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
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
function (t, N, pi0 = 0.5, prior.cauchy = TRUE, null.d = 0.2, 
    MCMC = FALSE, scale = 1) 
{
    nu = ifelse(prior.cauchy, 1, Inf)
    pi2 = pt(abs(null.d)/scale, nu) - pt(-abs(null.d)/scale, 
        nu)
    B1 = t.value.bf(t, N, prior.cauchy = prior.cauchy)
    B2 = area.bf(t, N, df.prior = nu, null.d = null.d, MCMC = MCMC)
    B3 = pi0 * pi2 * B1 * B2 + pi0 * B1 - pi0 * pi2 * B1 + B2 - 
        pi0 * B2
    B3
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
