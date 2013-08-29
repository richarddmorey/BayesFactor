##' Using the classical t test statistic for a one- or two-sample design, this 
##' function computes the corresponding Bayes factor test.
##' 
##' This function can be used to compute the Bayes factor corresponding to a 
##' one-sample, a paired-sample, or an independent-groups t test, using the 
##' classical t statistic.  It can be used when you don't have access to the 
##' full data set for analysis by \code{\link{ttestBF}}, but you do have the 
##' test statistic.
##' 
##' For details about the model, see the help for \code{\link{ttestBF}}, and the
##' references therein.
##' 
##' The Bayes factor is computed via Gaussian quadrature.
##' @title Use t statistic to compute Bayes factor for one- and two- sample designs
##' @param t classical t statistic
##' @param n1 size of first group (or only group, for one-sample tests)
##' @param n2 size of second group, for independent-groups tests
##' @param nullInterval optional vector of length 2 containing lower and upper bounds of an interval hypothesis to test, in standardized units
##' @param rscale numeric prior scale
##' @return If \code{nullInterval} is defined, then two Bayes factors will be
##'   computed: The Bayes factor for the interval against the null hypothesis 
##'   that the standardized effect is 0, and the corresponding Bayes factor for 
##'   the compliment of the interval. For each Bayes factor, a vector of length
##'   2 containing the computed log(e) Bayes factor (against the point null),
##'   along with a proportional error estimate on the Bayes factor is returned.
##' @author Richard D. Morey (\email{richarddmorey@@gmail.com}) and Jeffrey N. 
##'   Rouder (\email{rouderj@@missouri.edu})
##' @keywords htest
##' @export
##' @references Morey, R. D. & Rouder, J. N. (2011). Bayes Factor Approaches for
##'   Testing Interval Null Hypotheses. Psychological Methods, 16, 406-419
##'   
##'   Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., & Iverson, G. 
##'   (2009). Bayesian t-tests for accepting and rejecting the null hypothesis. 
##'   Psychonomic Bulletin & Review, 16, 752-760
##' @seealso \code{\link{integrate}}, \code{\link{t.test}}; see 
##'   \code{\link{ttestBF}} for the intended interface to this function, using 
##'   the full data set.
##' @examples 
##' ## Classical example: Student's sleep data
##' data(sleep)
##' plot(extra ~ group, data = sleep)
##' 
##' ## t.test() gives a t value of -4.0621
##' t.test(extra ~ group, data = sleep, paired=TRUE)
##' ## Gives a Bayes factor of about 15
##' ## in favor of the alternative hypothesis
##' result <- ttest.tstat(t = -4.0621, n1 = 10)
##' exp(result[['bf']])

ttest.tstat=function(t,n1,n2=0,nullInterval=NULL,rscale="medium")
{
  if(n2){
    rscale = rpriorValues("ttestTwo",,rscale)
  }else{
    rscale = rpriorValues("ttestOne",,rscale)
  }
  
  nu=ifelse(n2==0 | is.null(n2),n1-1,n1+n2-2)
  n=ifelse(n2==0 | is.null(n2),n1,(n1*n2)/(n1+n2))
  r2=rscale^2
  log.marg.like.0= -(nu+1)/2 * log(1+t^2/(nu))
  if(is.null(nullInterval)){
    integral = integrate(t.joint,lower=0,upper=Inf,t=t,n=n,nu=nu,r2=r2)  
    marg.like.1 = integral$value
    prop.error = integral$abs.error / marg.like.1
    lbf = log(marg.like.1) - log.marg.like.0
  }else{
    areabf = ttestAreaNull(t, n1, n2, nullInterval=nullInterval, rscale=rscale)
    lbf = areabf$bf
    prop.error = areabf$properror
  }
  return(list(bf = lbf, properror=prop.error))
}
