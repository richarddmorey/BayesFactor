\name{area.bf.integrate}
\Rdversion{1.1}
\alias{area.bf.integrate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
area.bf.integrate(t, N, df.prior = 1, null.d = 0.2, MCMC = FALSE, eta2.lo = 0, eta2.up = Inf, d.lo = -Inf, d.up = Inf, r.scale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
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
  \item{eta2.lo}{
%%     ~~Describe \code{eta2.lo} here~~
}
  \item{eta2.up}{
%%     ~~Describe \code{eta2.up} here~~
}
  \item{d.lo}{
%%     ~~Describe \code{d.lo} here~~
}
  \item{d.up}{
%%     ~~Describe \code{d.up} here~~
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
function (t, N, df.prior = 1, null.d = 0.2, MCMC = FALSE, eta2.lo = 0, 
    eta2.up = Inf, d.lo = -Inf, d.up = Inf, r.scale = 1) 
{
    null.d = abs(null.d)
    priorarea = pt(null.d/r.scale, df.prior) - pt(-null.d/r.scale, 
        df.prior)
    priorodds = priorarea/(1 - priorarea)
    postarea.null = integrate(vect.margd.post, -null.d, null.d, 
        t = t, N = N, df.prior = df.prior, r.scale = r.scale, 
        rel.tol = max(50 * .Machine$double.eps, 5e-29))$value
    postarea.total = integrate(vect.margd.post, d.lo, d.up, t = t, 
        N = N, df.prior = df.prior, r.scale = r.scale, rel.tol = max(50 * 
            .Machine$double.eps, 5e-29))$value
    post.d = function(x) vect.margd.post(x, t, N, df.prior = 1, 
        eta2.lo = eta2.lo, eta2.up = eta2.up, r.scale = r.scale)/postarea.total
    BF = (postarea.null/(postarea.total - postarea.null))/priorodds
    output = list(bf.area = BF, bfMCMC95 = 0, post.d = post.d, 
        t = t, N = N, nu = df.prior, iterations = 0, acc.rate.delta = 0, 
        acc.rate.eta2 = 0, chains = 0)
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
