##' Using the classical t test statistic for multiple  one- or two-sample designed-experiments, this 
##' function computes the corresponding meta-analytic Bayes factor test.
##' 
##' This function can be used to compute a meta-analytic Bayes factor corresponding to a 
##' one-sample, a paired-sample, or an independent-groups t test, using the 
##' classical t statistic. 
##' 
##' For details about the model, see the references below, and the help for \code{\link{ttestBF}}, and the
##' references therein.
##' 
##' The Bayes factor is computed via Gaussian quadrature.
##' @title Use t statistic to compute Bayes factor for one- and two- sample designs
##' @param t vector of classical t statistics
##' @param n1 vector of sizes of first group (or only group, for one-sample tests)
##' @param n2 vector of size of second group, for independent-groups tests
##' @param nullInterval optional vector of length 2 containing lower and upper bounds of an interval hypothesis to test, in standardized units
##' @param rscale numeric prior scale
##' @return If \code{nullInterval} is defined, then two Bayes factors will be
##'   computed: The Bayes factor for the interval against the null hypothesis 
##'   that the standardized effect is 0, and the corresponding Bayes factor for 
##'   the compliment of the interval. For each Bayes factor, a vector of length
##'   2 containing the computed log(e) Bayes factor (against the point null),
##'   along with a proportional error estimate on the Bayes factor is returned.
##' @author Richard D. Morey (\email{richarddmorey@@gmail.com})
##' @keywords htest
##' @export
##' @references Morey, R. D. & Rouder, J. N. (2011). Bayes Factor Approaches for
##'   Testing Interval Null Hypotheses. Psychological Methods, 16, 406-419
##'   
##'   Rouder, J. N. & Morey, R. D. (2011). A Bayes Factor Meta-Analysis of Bem's ESP Claim.
##'   Psychonomic Bulletin & Review, 18, 682-689  
##' @seealso \code{\link{integrate}}, \code{\link{t.test}}; see 
##'   \code{\link{ttest.tstat}} for the single-study version of this function
##' @examples 
##' ## Bem's (2010) data (see Rouder & Morey, 2011)
##' t=c(-.15,2.39,2.42,2.43)
##' N=c(100,150,97,99)
##' 
##' ## Using rscale=1 to be consistent with Rouder & Morey (2011)
##' exp(meta.ttest.tstat(t,n1=N,rscale=1)$bf) 

meta.ttest.tstat <- function(t,n1,n2=NULL,nullInterval=NULL,rscale="medium")
{
  
  if( ((length(n1) != length(n2)) & !is.null(n2)) | (length(n1) != length(t))){
    stop("Number of t statistics must equal number of sample sizes.")  
  }
  if(!is.null(n2)){
    rscale = BayesFactor:::rpriorValues("ttestTwo",,rscale)
  }else{
    rscale = BayesFactor:::rpriorValues("ttestOne",,rscale)
  }
  if(is.null(n2)){
    nu = n1 - 1  
    n <- n1
  }else{
    nu = n1 + n2 - 2
    n <- n1 * n2 / (n1 + n2)
  }
  
  log.marg.like.0 = sum(dt(t, df = nu, log=TRUE))
  if(is.null(nullInterval)){
    integral = suppressWarnings(integrate(t.joint.meta,lower=-Inf,upper=Inf,t=t,n=n,nu=nu,rscale=rscale,const=log.marg.like.0))  
    marg.like.1 = integral$value
    prop.error = integral$abs.error / marg.like.1
    lbf = log(marg.like.1)
  }else{
    areabf = ttestAreaNull.meta(t, n1, n2, nullInterval=nullInterval, rscale=rscale)
    lbf = areabf$bf
    prop.error = areabf$properror
  }
  return(list(bf = lbf, properror=prop.error))  
}

t.joint.meta = Vectorize(function(delta,priorFunc,t,n,nu,rscale, const=0)
  exp(sum(dt(t, df = nu, ncp = delta * sqrt(n), log = TRUE)) + 
        dcauchy(delta, scale=rscale, log=TRUE) - const),
  "delta")

ttestAreaNull.meta <- function(t, n1, n2=0, nullInterval=c(-.2,.2), rscale, safeInt = .9999)
{
  if(is.null(n2)){
    nu = n1 - 1  
    N <- n1
  }else{
    nu = n1 + n2 - 2
    N <- n1 * n2 / (n1 + n2)
  }
  
  nullInterval = range(nullInterval)
  
  priorOdds = diff(pcauchy(nullInterval,scale=rscale))
   
  
  ## Use normal approximation to guess good integration limits
  meanDelta = sum(N*t/sqrt(N))/sum(N)
  sdDelta = 1/sqrt(sum(N))
  safeRange = meanDelta + c(-1,1) * qt(1-(1-safeInt)/2,df=sum(nu))*sdDelta
  
  priorOdds = diff(pcauchy(nullInterval,scale=rscale))
  
  nullInterval[1] = max(nullInterval[1],safeRange[1]) 
  nullInterval[2] = min(nullInterval[2],safeRange[2]) 
    
  allIntegral = suppressWarnings(integrate(Vectorize(function(delta,tstat,n1,nu,rscale){
    exp(sum(dt(tstat,df = nu, ncp = delta * sqrt(n1), log=TRUE)) + dcauchy(delta, scale=rscale, log=TRUE))
  },"delta"), -Inf, Inf, tstat=t,n1=N,nu=nu, rscale=rscale)[[1]])
  
  areaIntegral = suppressWarnings(integrate(Vectorize(function(delta,tstat,n1,nu,rscale,const=1){
    exp(sum(dt(tstat,df = nu, ncp = delta * sqrt(n1), log=TRUE)) + dcauchy(delta, scale=rscale, log=TRUE) - log(const))
  },"delta"), nullInterval[1], nullInterval[2], tstat=t,n1=N,nu=nu,rscale=rscale,const=allIntegral))
  
  
  # encompassing vs point null
  vsNull = meta.ttest.tstat(t, n1, n2, rscale=rscale)
  
  val = areaIntegral[[1]]
  err = areaIntegral[[2]]
  
  err = err / val 
  err = sqrt(err^2 + vsNull[['properror']]^2)
  val = log(val) -  log(priorOdds) + vsNull[['bf']]
  
  complArea = 1-areaIntegral[[1]]
  errCompl = areaIntegral[[2]] / complArea
  errCompl = sqrt(errCompl^2 + vsNull[['properror']]^2)  
  valCompl = log(complArea) - log(1-priorOdds) + vsNull[['bf']]
  
  return(
    list(
      bf = c(null=val,alt=valCompl),
      properror = c(null=err,alt=errCompl)
    ))
}




