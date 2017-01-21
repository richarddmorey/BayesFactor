##' Bayes factors or posterior samples for binomial, geometric, or neg. binomial data.
##'
##' Given count data modeled as a binomial, geometric, or negative binomial random variable,
##' the Bayes factor provided by \code{proportionBF} tests the null hypothesis that
##' the probability of a success is \eqn{p_0}{p_0} (argument \code{p}). Specifically,
##' the Bayes factor compares two hypotheses: that the probability is \eqn{p_0}{p_0}, or
##' probability is not \eqn{p_0}{p_0}. Currently, the default alternative is that
##' \deqn{\lambda~logistic(\lambda_0,r)} where
##' \eqn{\lambda_0=logit(p_0)}{lambda_0=logit(p_0)} and
##' \eqn{\lambda=logit(p)}{lambda=logit(p)}. \eqn{r}{r} serves as a prior scale parameter.
##'
##' For the \code{rscale} argument, several named values are recognized:
##' "medium", "wide", and "ultrawide". These correspond
##' to \eqn{r} scale values of \eqn{1/2}{1/2}, \eqn{\sqrt{2}/2}{sqrt(2)/2}, and 1,
##' respectively.
##'
##' The Bayes factor is computed via Gaussian quadrature, and posterior
##' samples are drawn via independence Metropolis-Hastings.
##' @title Function for Bayesian analysis of proportions
##' @param y a vector of successes
##' @param N a vector of total number of observations
##' @param p the null value for the probability of a success to be tested against
##' @param rscale prior scale.  A number of preset values can be given as
##'   strings; see Details.
##' @param nullInterval optional vector of length 2 containing
##' lower and upper bounds of an interval hypothesis to test, in probability units
##' @param posterior if \code{TRUE}, return samples from the posterior instead
##'   of Bayes factor
##' @param callback callback function for third-party interfaces
##' @param ... further arguments to be passed to or from methods.
##' @return If \code{posterior} is \code{FALSE}, an object of class
##'   \code{BFBayesFactor} containing the computed model comparisons is
##'   returned. If \code{nullInterval} is defined, then two Bayes factors will
##'   be computed: The Bayes factor for the interval against the null hypothesis
##'   that the probability is \eqn{p_0}{p0}, and the corresponding Bayes factor for
##'   the compliment of the interval.
##'
##'   If \code{posterior} is \code{TRUE}, an object of class \code{BFmcmc},
##'   containing MCMC samples from the posterior is returned.
##' @export
##' @keywords htest
##' @author Richard D. Morey (\email{richarddmorey@@gmail.com})
##' @examples
##' bf = proportionBF(y = 15, N = 25, p = .5)
##' bf
##' ## Sample from the corresponding posterior distribution
##' samples =proportionBF(y = 15, N = 25, p = .5, posterior = TRUE, iterations = 10000)
##' plot(samples[,"p"])
##' @seealso \code{\link{prop.test}}

proportionBF <- function(y, N, p, rscale = "medium", nullInterval = NULL, posterior=FALSE, callback = function(...) as.integer(0), ...)
{
  if (p >= 1 || p <= 0)
    stop('p must be between 0 and 1', call.=FALSE)

  if(!is.null(nullInterval)){
    if(any(nullInterval<0) | any(nullInterval>1)) stop("nullInterval endpoints must be in [0,1].")
    nullInterval = range(nullInterval)
  }
  rscale = rpriorValues("proptest",,rscale)

  if( length(p) > 1 ) stop("Only a single null allowed (length(p) > 1).")
  if( length(y) != length(N) ) stop("Length of y and N must be the same.")
  if( any(y>N) | any(y < 0) ) stop("Invalid data (y>N or y<0).")
  if( any( c(y,N)%%1 != 0 ) ) stop("y and N must be integers.")

  hypNames = makePropHypothesisNames(rscale, nullInterval, p)

  mod1 = BFproportion(type = "logistic",
                 identifier = list(formula = "p =/= p0", nullInterval = nullInterval, p0 = p),
                 prior=list(rscale=rscale, nullInterval = nullInterval, p0 = p),
                 shortName = hypNames$shortName,
                 longName = hypNames$longName
  )

  data = data.frame(y = y, N = N)

  checkCallback(callback,as.integer(0))

  if(posterior)
    return(posterior(mod1, data = data, callback = callback, ...))

  bf1 = compare(numerator = mod1, data = data)

  if(!is.null(nullInterval)){
    mod2 = mod1
    attr(mod2@identifier$nullInterval, "complement") = TRUE
    attr(mod2@prior$nullInterval, "complement") = TRUE
    hypNames = makePropHypothesisNames(rscale, mod2@identifier$nullInterval,p)
    mod2@shortName = hypNames$shortName
    mod2@longName = hypNames$longName

    bf2 = compare(numerator = mod2, data = data)
    checkCallback(callback,as.integer(1000))
    return(c(bf1, bf2))
  }else{
    checkCallback(callback,as.integer(1000))
    return(c(bf1))
  }
}
