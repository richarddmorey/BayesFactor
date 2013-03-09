##' This function computes Bayes factors, or samples from the posterior, for 
##' one- and two-sample designs.
##' 
##' The Bayes factor provided by \code{ttestBF} tests the null hypothesis that 
##' the mean (or mean difference) of a normal population is \eqn{\mu_0}{mu0} 
##' (argument \code{mu}). Specifically, the Bayes factor compares two 
##' hypotheses: that the standardized effect size is 0, or that the standardized
##' effect size is not 0. For one-sample tests, the standardized effect size is 
##' \eqn{(\mu-\mu_0)/\sigma}{(mu-mu0)/sigma}; for two sample tests, the 
##' standardized effect size is \eqn{(\mu_2-\mu_1)/\sigma}{(mu2-mu1)/sigma}.
##' 
##' A noninformative Jeffreys prior is placed on the variance of the normal 
##' population, while a Cauchy prior is placed on the standardized effect size. 
##' The \code{rscale} argument controls the scale of the prior distribution, 
##' with \code{rscale=1} yielding a standard Cauchy prior. See the references 
##' below for more details.
##' 
##' For the \code{rscale} argument, several named values are recognized: 
##' "medium" corresponds to \eqn{r=\sqrt{2}/2}{r=sqrt(2)/2}; "wide" corresponds 
##' to \eqn{r=1}{r=1}.
##' 
##' The Bayes factor is computed via Gaussian quadrature.
##' @title Function for Bayesian analysis of one- and two-sample designs
##' @param x a vector of observations for the first (or only) group
##' @param y a vector of observations for the second group (or condition, for 
##'   paired)
##' @param formula for independent-group designs, a (optional) formula 
##'   describing the model
##' @param mu for one-sample and paired designs, the null value of the mean (or 
##'   mean difference)
##' @param nullInterval optional vector of length 2 containing lower and upper bounds of an interval hypothesis to test, in standardized units
##' @param paired if \code{TRUE}, observations are paired
##' @param data for use with \code{formula}, a data frame containing all the 
##'   data
##' @param rscale prior scale.  A number of preset values can be given as 
##'   strings; see Details.
##' @param posterior if \code{TRUE}, return samples from the posterior instead 
##'   of Bayes factor
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
##' @references Morey, R. D., Rouder, J. N., Pratte, M. S., & Speckman, P. L. 
##'   (2011). Using MCMC chain outputs to efficiently estimate Bayes factors. 
##'   Journal of Mathematical Psychology, 55, 368-378
##'   
##'   Morey, R. D. \& Rouder, J. N. (2011). Bayes Factor Approaches for Testing 
##'   Interval Null Hypotheses. Psychological Methods, 16, 406-419
##'   
##'   Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., & Iverson, G. 
##'   (2009). Bayesian t-tests for accepting and rejecting the null hypothesis. 
##'   Psychonomic Bulletin & Review, 16, 752-760
##'   
##'   Perception and Cognition Lab (University of Missouri): Bayes factor 
##'   calculators. \url{http://pcl.missouri.edu/bayesfactor}
##' @note The default priors have scale has changed from 1 to \eqn{\sqrt{2}/2} for the 
##'   two-sample t test, and 1/2 for the one-sample t test. The 
##'   factor of \eqn{\sqrt{2}} in the two-sample t test is to be consistent 
##'   with Morey et al. (2011) and 
##'   Rouder et al. (2012), and the factor of \eqn{1/2} in both is to better scale the 
##'   expected effect sizes; the previous scaling put more weight on larger 
##'   effect sizes. To obtain the same Bayes factors as Rouder et al. (2009), 
##'   change the prior scale to 1.
##' @examples
##' ## Sleep data from t test example
##' data(sleep)
##' plot(extra ~ group, data = sleep)
##' 
##' ## paired t test
##' ttestBF(x = sleep$extra[sleep$group==1], y = sleep$extra[sleep$group==2], paired=TRUE)
##' 
##' ## Sample from the corresponding posterior distribution
##' samples = ttestBF(x = sleep$extra[sleep$group==1], y = sleep$extra[sleep$group==2], paired=TRUE, posterior = TRUE, iterations = 1000)
##' plot(samples[,"mu"])
##' @seealso \code{\link{integrate}}, \code{\link{t.test}}

ttestBF <- function(x, y = NULL, formula = NULL, mu = 0, nullInterval = NULL, 
                    paired = FALSE, data = NULL, rscale="medium", posterior=FALSE, ...){
  
  
  if( (is.null(formula) & is.null(y)) | (!is.null(y) & paired) ){
    if(paired){
      # check that the two vectors have same length
      if(length(x)!=length(y)) stop("Length of x and y must be the same if paired=TRUE.")
      x = x - y
    }
    rscale = rpriorValues("ttestOne",,rscale)
    if(is.null(nullInterval)){
      modFull = BFoneSample(type = "JZS", 
                            identifier = list(formula = "y ~ 1"), 
                            prior=list(rscale=rscale, mu=mu),
                            shortName = paste("Alt., r=",round(rscale,3),sep=""),
                            longName = paste("Alternative, r = ",rscale,", mu =/= ",mu, sep="")
      )
      
      if(posterior){
        chains = posterior(modFull,data = data.frame(y=x), ...)
        return(chains)
      }else{
        bf = compare(numerator = modFull, data = data.frame(y=x))
        return(bf)
      }
    }else{
      nullInterval = range(nullInterval)
      modInterval = BFoneSample(type = "JZS", 
                                 identifier = list(formula = "y ~ 1",nullInterval = nullInterval), 
                                 prior=list(rscale=rscale, mu=mu, nullInterval = nullInterval),
                                 shortName = paste("Alt., r=",round(rscale,3)," ",nullInterval[1],"<d<",nullInterval[2],sep=""),
                                 longName = paste("Alternative, r = ",rscale,", mu =/= ",mu, " ",nullInterval[1],"<d<",nullInterval[2],sep="")
      )      
      if(posterior){
        chains = posterior(modInterval, data = data.frame(y=x), ...)
        return(chains)
      }else{
        bf = compare(numerator = modInterval, data = data.frame(y=x))
        return(bf)
      }
    }
  }else if(!is.null(y) & !paired){
    data = data.frame(y = c(x,y), 
                      group = factor(c(rep("x",length(x)),rep("y",length(y))))
                      )
    formula = y ~ group
  }
  if(!is.null(formula)){ # formula
    if(paired) stop("Cannot use 'paired' with formula.")
    if(is.null(data)) stop("'data' needed for formula.")
    
    ivs = attr(terms(formula, data = data),"term.labels")
    if(length(ivs) > 1) stop("Only one independent variable allowed for t test.")
    dataTypes = "fixed"
    names(dataTypes) = ivs
    if(mu != 0) stop("Use of nonzero null hypothesis not implemented for independent samples test.")
    
    rscale = rpriorValues("ttestTwo",,rscale)
    
    if(is.null(nullInterval)){
      numerator = BFindepSample(type = "JZS", 
                              identifier = list(formula = stringFromFormula(formula)), 
                              prior=list(rscale=rscale, mu=mu),
                              shortName = paste("Alt., r=",round(rscale,3),sep=""),
                              longName = paste("Alternative, r = ",rscale,", mu =/= ",mu,sep="")
      )
    }else{
      nullInterval = range(nullInterval)
      numerator = BFindepSample(type = "JZS", 
                                identifier = list(formula = stringFromFormula(formula),nullInterval = nullInterval), 
                                prior=list(rscale=rscale, mu=mu, nullInterval = nullInterval),
                                shortName = paste("Alt., r=",round(rscale,3)," ",nullInterval[1],"<d<",nullInterval[2],sep=""),
                                longName = paste("Alternative, r = ",rscale,", mu =/= ",mu, " ",nullInterval[1],"<d<",nullInterval[2],sep="")
      )
    }

    if(posterior){
      chains = posterior(numerator, data = data, ...)
      return(chains)
    }else{
      bf = compare(numerator = numerator, data = data)
      return(bf)
    }
  }
}



