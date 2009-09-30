\name{PCLBF}
\Rdversion{1.1}
\alias{PCLBF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PCLBF(t, N1, N2 = 0, prior.cauchy = TRUE, null.d = 0.2, pi0 = 0.5, scale = 1, MCMC = FALSE, CIconf = 0.95, ...)
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
  \item{prior.cauchy}{
%%     ~~Describe \code{prior.cauchy} here~~
}
  \item{null.d}{
%%     ~~Describe \code{null.d} here~~
}
  \item{pi0}{
%%     ~~Describe \code{pi0} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{MCMC}{
%%     ~~Describe \code{MCMC} here~~
}
  \item{CIconf}{
%%     ~~Describe \code{CIconf} here~~
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
function (t, N1, N2 = 0, prior.cauchy = TRUE, null.d = 0.2, pi0 = 0.5, 
    scale = 1, MCMC = FALSE, CIconf = 0.95, ...) 
{
    BF = new("PCLBF")
    BF@t.value = t
    BF@N1 = N1
    BF@N2 = N2
    BF@prior.cauchy = prior.cauchy
    BF@MCMC = MCMC
    BF@pi0 = pi0
    BF@null.d = null.d
    BF@prior.scale = scale
    t.df = ifelse(N2 == 0 | is.null(N2), N1 - 1, N1 + N2 - 2)
    my.N = ifelse(N2 == 0 | is.null(N2), N1, (N1 * N2)/(N1 + 
        N2))
    BF@num.samples = ifelse(N2 == 0 | is.null(N2), 1, 2)
    nu = ifelse(prior.cauchy, 1, Inf)
    pi2 = pt(abs(null.d)/scale, nu) - pt(-abs(null.d)/scale, 
        nu)
    BF@JZS = t.value.bf(t, N1, N2, prior.cauchy = prior.cauchy, 
        sd = scale)
    area = area.bf(t, N1, N2, nu, null.d, MCMC = MCMC, r.scale = scale, 
        ...)
    BF@NOH = area$bf.area
    BF@marginal.d = area$post.d
    BF@MCMC.iters = area$iterations
    BF@MCMC.95CI = area$bfMCMC95
    BF@MCMC.acceptance = c(area$acc.rate.delta, area$acc.rate.eta2)
    BF@MCMC.chains = area$chains
    BF@freq.CI = freqCIdelta(t, N1, alpha = 1 - CIconf)
    BF@freq.conf = CIconf
    BF@freq.p = 2 * (1 - pt(abs(t), t.df))
    BF@hybrid = pi0 * pi2 * BF@JZS * BF@NOH + pi0 * BF@JZS - 
        pi0 * pi2 * BF@JZS + BF@NOH - pi0 * BF@NOH
    return(BF)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
