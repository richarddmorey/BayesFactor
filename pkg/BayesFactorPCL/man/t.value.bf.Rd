\name{t.value.bf}
\Rdversion{1.1}
\alias{t.value.bf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
t.value.bf(t, n1, n2 = 0, sd = 1, prior.cauchy = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{n1}{
%%     ~~Describe \code{n1} here~~
}
  \item{n2}{
%%     ~~Describe \code{n2} here~~
}
  \item{sd}{
%%     ~~Describe \code{sd} here~~
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
function (t, n1, n2 = 0, sd = 1, prior.cauchy = T) 
{
    nu = ifelse(n2 == 0 | is.null(n2), n1 - 1, n1 + n2 - 2)
    n = ifelse(n2 == 0 | is.null(n2), n1, (n1 * n2)/(n1 + n2))
    r2 = sd * sd
    marg.like.0 = (1 + t^2/(nu))^(-(nu + 1)/2)
    marg.like.1 = ifelse(prior.cauchy, integrate(t.joint, lower = 0, 
        upper = Inf, t = t, n = n, nu = nu, r2 = r2)$value, (1 + 
        n * r2)^(-0.5) * (1 + t^2/((1 + n * r2) * (nu)))^(-(nu + 
        1)/2))
    return(marg.like.0/marg.like.1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
