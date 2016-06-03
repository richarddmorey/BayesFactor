##' Bayes factors or posterior samples for correlations.
##' 
##' 
##' For the \code{rscale} argument, several named values are recognized: 
##' "medium.narrow", "medium", "wide", and "ultrawide". These correspond 
##' to \eqn{r} scale values of \eqn{1/\sqrt(27)}{1/sqrt(27)}, \eqn{1/3){1/3}, 
##' \eqn{1/\sqrt(3)}{1/sqrt(3)} and 1, respectively.
##' 
##' The Bayes factor is computed via several different methods.
##' @title Function for Bayesian analysis of correlations
##' @param y first continuous variable
##' @param x second continuous variable
##' @param rscale prior scale.  A number of preset values can be given as 
##'   strings; see Details.
##' @param nullInterval optional vector of length 2 containing 
##' lower and upper bounds of an interval hypothesis to test, in correlation units
##' @param posterior if \code{TRUE}, return samples from the posterior instead 
##'   of Bayes factor
##' @param callback callback function for third-party interfaces 
##' @param ... further arguments to be passed to or from methods.
##' @return If \code{posterior} is \code{FALSE}, an object of class 
##'   \code{BFBayesFactor} containing the computed model comparisons is 
##'   returned. If \code{nullInterval} is defined, then two Bayes factors will
##'   be computed: The Bayes factor for the interval against the null hypothesis
##'   that the probability is 0, and the corresponding Bayes factor for
##'   the complement of the interval.
##'   
##'   If \code{posterior} is \code{TRUE}, an object of class \code{BFmcmc},
##'   containing MCMC samples from the posterior is returned.
##' @export
##' @keywords htest
##' @author Richard D. Morey (\email{richarddmorey@@gmail.com})
##' @examples 
##' bf = correlationBF(y = iris$Sepal.Length, x = iris$Sepal.Width)
##' bf
##' ## Sample from the corresponding posterior distribution
##' samples = correlationBF(y = iris$Sepal.Length, x = iris$Sepal.Width, posterior = TRUE, iterations = 10000)
##' plot(samples[,"rho"])
##' @seealso \code{\link{cor.test}}
##' @references Ly, A., Verhagen, A. J. & Wagenmakers, E.-J. (2015). 
##' Harold Jeffreys's Default Bayes Factor Hypothesis Tests: Explanation, Extension, and Application in Psychology.
##' Journal of Mathematical Psychology, Available online 28 August 2015, http://dx.doi.org/10.1016/j.jmp.2015.06.004.

correlationBF <- function(y, x, rscale = "medium", nullInterval = NULL, posterior=FALSE, callback = function(...) as.integer(0), ...)
{
  if(!is.null(nullInterval)){
    if(any(nullInterval< -1) | any(nullInterval>1)) stop("nullInterval endpoints must be in [-1,1].")
    nullInterval = range(nullInterval)
  }
  rscale = rpriorValues("correlation",,rscale)
  
  if( length(y) != length(x) ) stop("Length of y and x must be the same.")
  if(!is.null(nullInterval))
    if(any(abs(nullInterval)>1))
      stop("Invalid interval hypothesis; endpoints must be in [-1,1].")
  if( length(y)<3 )  stop("N must be >2.")
  n = length(x) - sum(is.na(y) | is.na(x))
  if(n < 3) stop("Need at least 3 complete observations.")
  
  hypNames = makeCorrHypothesisNames(rscale, nullInterval)
  
  mod1 = BFcorrelation(type = "Jeffreys-beta*", 
                 identifier = list(formula = "rho =/= 0", nullInterval = nullInterval), 
                 prior=list(rscale=rscale, nullInterval = nullInterval),
                 shortName = hypNames$shortName,
                 longName = hypNames$longName
  )      
  
  data = data.frame(y = y, x = x)
  
  checkCallback(callback,as.integer(0))
  
  if(posterior)
    return(posterior(mod1, data = data, callback = callback, ...))
  
  bf1 = compare(numerator = mod1, data = data)
  
  if(!is.null(nullInterval)){
    mod2 = mod1
    attr(mod2@identifier$nullInterval, "complement") = TRUE
    attr(mod2@prior$nullInterval, "complement") = TRUE
    hypNames = makeCorrHypothesisNames(rscale, mod2@identifier$nullInterval)
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
