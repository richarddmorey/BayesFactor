\name{area.bf.mcmc}
\Rdversion{1.1}
\alias{area.bf.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
area.bf.mcmc(t, N, nu = 1, null.d = 0.2, iters = 10000, metscale = 2, pnull = 0.5, chains = FALSE, r.scale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{nu}{
%%     ~~Describe \code{nu} here~~
}
  \item{null.d}{
%%     ~~Describe \code{null.d} here~~
}
  \item{iters}{
%%     ~~Describe \code{iters} here~~
}
  \item{metscale}{
%%     ~~Describe \code{metscale} here~~
}
  \item{pnull}{
%%     ~~Describe \code{pnull} here~~
}
  \item{chains}{
%%     ~~Describe \code{chains} here~~
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
function (t, N, nu = 1, null.d = 0.2, iters = 10000, metscale = 2, 
    pnull = 0.5, chains = FALSE, r.scale = 1) 
{
    null.d = abs(null.d)
    pnull.cor = ifelse(nu == Inf, pnorm(null.d/r.scale) - pnorm(-null.d/r.scale), 
        pt(null.d/r.scale, nu) - pt(-null.d/r.scale, nu))
    area.cor = pnull.cor/(1 - pnull.cor)
    prior.cor = pnull/(1 - pnull)
    a = N/2
    b = (N - 1) * 1/2
    sdfullcond.eta2 = b/((a - 1) * sqrt(a - 2))
    sdfullcond.delta = 1/sqrt(N)
    metscale.eta2 = sdfullcond.eta2 * metscale
    metscale.delta = sdfullcond.delta * metscale
    eta2.samp = 1:iters * NA
    delta.samp = 1:iters * NA
    eta2.samp[1] = 1
    delta.samp[1] = t/sqrt(N)
    for (i in 2:iters) {
        eta2.samp[i] = samp.eta2.t(eta2.samp[i - 1], delta.samp[i - 
            1], t, N, metscale.eta2)
        delta.samp[i] = samp.delta.t(delta.samp[i - 1], eta2.samp[i], 
            t, N, metscale.delta, nu = nu, r = r.scale)
    }
    postmean.delta = mean(delta.samp)
    acc.rate.delta = mean(diff(delta.samp) != 0)
    acc.rate.eta2 = mean(diff(eta2.samp) != 0)
    postprobnull = mean(delta.samp > -null.d & delta.samp < null.d)
    postoddsnull = postprobnull/(1 - postprobnull)
    bf.area = (postoddsnull/area.cor)/prior.cor
    innull = delta.samp > -null.d & delta.samp < null.d
    phat = mean(innull)
    mcmc95 = phat + c(-1, 1) * 1.96 * sqrt(phat * (1 - phat))/sqrt(effectiveSize(as.numeric(innull)))
    bfCI95 = ((mcmc95/(1 - mcmc95))/area.cor)/prior.cor
    post.d = approxfun(density(delta.samp)$x, density(delta.samp)$y)
    if (chains) {
        my.chains = as.mcmc(data.frame(delta.samp, eta2.samp, 
            as.numeric(innull)))
    }
    else {
        my.chains = 0
    }
    output = list(bf.area = bf.area, bfMCMC95 = bfCI95, post.d = post.d, 
        t = t, N = N, nu = nu, iterations = iters, acc.rate.delta = acc.rate.delta, 
        acc.rate.eta2 = acc.rate.eta2, chains = my.chains)
    output
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
