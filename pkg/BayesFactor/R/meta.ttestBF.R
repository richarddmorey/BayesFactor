##' This function computes mata-analytic Bayes factors, or samples from the posterior, for 
##' one- and two-sample designs where multiple t values have been observed.
##' 
##' The Bayes factor provided by \code{meta.ttestBF} tests the null hypothesis that 
##' the true effect size (or alternatively, the noncentrality parameters) underlying a 
##' set of t statistics is 0. Specifically, the Bayes factor compares two 
##' hypotheses: that the standardized effect size is 0, or that the standardized
##' effect size is not 0. Note that there is assumed to be a single, common effect size 
##' \eqn{\delta}{delta} underlying all t statistics. For one-sample tests, the standardized effect size is 
##' \eqn{(\mu-\mu_0)/\sigma}{(mu-mu0)/sigma}; for two sample tests, the 
##' standardized effect size is \eqn{(\mu_2-\mu_1)/\sigma}{(mu2-mu1)/sigma}.
##' 
##' A Cauchy prior is placed on the standardized effect size. 
##' The \code{rscale} argument controls the scale of the prior distribution, 
##' with \code{rscale=1} yielding a standard Cauchy prior. See the help for 
##' \code{\link{ttestBF}} and the references below for more details.
##' 
##' The Bayes factor is computed via Gaussian quadrature. Posterior samples are 
##' drawn via independent-candidate Metropolis-Hastings.
##' @title Function for Bayesian analysis of one- and two-sample designs
##' @param t a vector of t statistics
##' @param n1 a vector of sample sizes for the first (or only) condition
##' @param n2 a vector of sample sizes. If \code{NULL}, a one-sample design is assumed
##' @param nullInterval optional vector of length 2 containing lower and upper bounds of 
##' an interval hypothesis to test, in standardized units
##' @param rscale prior scale.  A number of preset values can be given as 
##'   strings; see Details.
##' @param posterior if \code{TRUE}, return samples from the posterior instead 
##'   of Bayes factor
##' @param callback callback function for third-party interfaces
##' @param ... further arguments to be passed to or from methods.
##' @return If \code{posterior} is \code{FALSE}, an object of class 
##'   \code{BFBayesFactor} containing the computed model comparisons is 
##'   returned. If \code{nullInterval} is defined, then two Bayes factors will
##'   be computed: The Bayes factor for the interval against the null hypothesis
##'   that the standardized effect is 0, and the corresponding Bayes factor for
##'   the compliment of the interval.
##'   
##'   If \code{posterior} is \code{TRUE}, an object of class \code{BFmcmc},
##'   containing MCMC samples from the posterior is returned.
##' @export
##' @keywords htest
##' @author Richard D. Morey (\email{richarddmorey@@gmail.com})
##' @references Morey, R. D. & Rouder, J. N. (2011). Bayes Factor Approaches for Testing 
##'   Interval Null Hypotheses. Psychological Methods, 16, 406-419
##'   
##'   Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., & Iverson, G. 
##'   (2009). Bayesian t-tests for accepting and rejecting the null hypothesis. 
##'   Psychonomic Bulletin & Review, 16, 752-760
##'   
##'   Rouder, J. N. & Morey, R. D. (2011). A Bayes Factor Meta-Analysis of Bem's ESP Claim.
##'   Psychonomic Bulletin & Review, 18, 682-689  
##' @note To obtain the same Bayes factors as Rouder and Morey (2011), 
##'   change the prior scale to 1.
##' @examples 
##' ## Bem's (2010) data (see Rouder & Morey, 2011)
##' t=c(-.15,2.39,2.42,2.43)
##' N=c(100,150,97,99)
##' 
##' ## Using rscale=1 and one-sided test to be 
##' ## consistent with Rouder & Morey (2011)
##' bf = meta.ttestBF(t, N, rscale=1, nullInterval=c(0, Inf)) 
##' bf[1]
##' 
##' ## plot posterior distribution of delta, assuming alternative
##' ## turn off progress bar for example
##' samples = posterior(bf[1], iterations = 1000, progress = FALSE)
##' ## Note that posterior() respects the nullInterval
##' plot(samples)
##' summary(samples)
##' @seealso \code{\link{ttestBF}}

meta.ttestBF <- function(t, n1, n2 = NULL, nullInterval = NULL, rscale="medium", 
                         posterior=FALSE, callback = function(...) as.integer(0), ...)
{
  rscale = rpriorValues(ifelse(is.null(n2),"ttestOne","ttestTwo"),,rscale)
  hypNames = makeMetaTtestHypothesisNames(rscale, nullInterval)
  data = data.frame(t=t,n1=n1)
  if(!is.null(n2)) data$n2 = n2
  if(!is.null(nullInterval)) nullInterval = range(nullInterval)
  
  mod1 = BFmetat(type = "JZS", 
                        identifier = list(formula = "delta =/= 0", nullInterval = nullInterval), 
                        prior=list(rscale=rscale, nullInterval = nullInterval),
                        shortName = hypNames$shortName,
                        longName = hypNames$longName
  )      
  
  if(posterior)
    return(posterior(mod1, data = data, callback = callback, ...))
  
  bf1 = compare(numerator = mod1, data = data)
    
  if(!is.null(nullInterval)){
    mod2 = mod1
    attr(mod2@identifier$nullInterval, "complement") = TRUE
    attr(mod2@prior$nullInterval, "complement") = TRUE
    hypNames = makeMetaTtestHypothesisNames(rscale, mod2@identifier$nullInterval)
    mod2@shortName = hypNames$shortName
    mod2@longName = hypNames$longName
    
    bf2 = compare(numerator = mod2, data = data)
    return(c(bf1, bf2))
  }else{
    return(c(bf1))
  }  
}

