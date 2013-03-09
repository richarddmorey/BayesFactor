##' This function computes Bayes factors, or samples from the posterior, of
##' specific linear models (either ANOVA or regression).
##' 
##' This function provides an interface for computing Bayes factors  for
##' specific linear models against the intercept-only null; other tests may be
##' obtained by computing two models and dividing their Bayes factors. Specifics
##' about the priors for regression models can be found in the help for
##' \code{\link{regressionBF}}; likewise, details for ANOVA models can be found
##' in the help for \code{\link{anovaBF}}.
##' 
##' Currently, the function does not allow for general linear models, containing
##' both continuous and categorical predcitors, but this support will be added
##' in the future.
##' @title Function to compute Bayes factors for specific linear models
##' @param formula a formula containing all factors to include in the analysis 
##'   (see Examples)
##' @param data a data frame containing data for all factors in the formula
##' @param whichRandom a character vector specifying which factors are random
##' @param rscaleFixed prior scale for standardized, reduced fixed effects. A 
##'   number of preset values can be given as strings; see Details.
##' @param rscaleRandom prior scale for standardized random effects
##' @param rscaleCont prior scale for standardized slopes
##' @param posterior if \code{TRUE}, return samples from the posterior
##'   distribution instead of the Bayes factor
##' @param ... further arguments to be passed to or from methods.
##' @return If \code{posterior} is \code{FALSE}, an object of class
##'   \code{BFBayesFactor}, containing the computed model comparisons is
##'   returned. Otherwise, an object of class \code{BFmcmc}, containing MCMC
##'   samples from the posterior is returned.
##' @author Richard D. Morey (\email{richarddmorey@@gmail.com})
##' @export
##' @keywords htest
##' @examples
##' ## Puzzles data; see ?puzzles for details
##' data(puzzles)
##' ## Bayes factor of full model against null
##' bfFull = lmBF(RT ~ shape + color + shape:color + ID, data = puzzles, whichRandom = "ID")
##' 
##' ## Bayes factor of main effects only against null
##' bfMain = lmBF(RT ~ shape + color + ID, data = puzzles, whichRandom = "ID")
##' 
##' ## Compare the main-effects only model to the full model
##' bfMain / bfFull
##' 
##' ## sample from the posterior of the full model
##' samples = lmBF(RT ~ shape + color + shape:color + ID, data = puzzles, whichRandom = "ID", posterior = TRUE, iterations = 1000)
##'  
##' ## Aother way to sample from the posterior of the full model
##' samples2 = posterior(bfFull, iterations = 1000)
##' @seealso  \code{\link{regressionBF}} and \code{anovaBF} for 
##' testing many regression or ANOVA models simultaneously.

lmBF <- function(formula, data, whichRandom = NULL, rscaleFixed="medium",
                 rscaleRandom=1, rscaleCont=1, posterior=FALSE, ...)
{    
  checkFormula(formula, data, analysis="lm")
  dataTypes <- createDataTypes(formula, whichRandom = whichRandom, data = data, analysis="lm")
  rscales = list(fixed=rpriorValues("allNways","fixed",rscaleFixed), 
                 random=rpriorValues("allNways","random",rscaleRandom), 
                 continuous=rpriorValues("regression",,rscaleCont))
  
  numerator = BFlinearModel(type = "JZS", 
                            identifier = list(formula = stringFromFormula(formula)), 
                            prior=list(rscale=rscales),
                            dataTypes = dataTypes,
                            shortName = paste(stringFromFormula(formula[[3]]),sep=""),
                            longName = paste(stringFromFormula(formula),sep="")
  )
  
  if(posterior){
    chains = posterior(numerator, data = data, ...)
    return(chains)
  }else{
    bf = compare(numerator = numerator, data = data,...)
    return(bf)    
  }
}
