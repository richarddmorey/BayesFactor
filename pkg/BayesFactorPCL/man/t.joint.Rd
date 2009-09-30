\name{t.joint}
\Rdversion{1.1}
\alias{t.joint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
t.joint(g, t, n, nu, r2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{nu}{
%%     ~~Describe \code{nu} here~~
}
  \item{r2}{
%%     ~~Describe \code{r2} here~~
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
function (g, t, n, nu, r2) 
{
    t1 = -0.5 * log(1 + n * g * r2)
    t2 = (-(nu + 1)/2) * log(1 + t^2/((1 + n * g * r2) * (nu)))
    return(dinvgamma(g, 0.5, 0.5) * exp(t1 + t2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
