##' This function computes Bayes factors corresponding to restrictions on a full model.
##' 
##' See the help for \code{\link{anovaBF}} and \code{\link{anovaBF}} or details.
##' 
##' Models, priors, and methods of computation are provided in Rouder et al. 
##' (2012) and Liang et al (2008).
##' 
##' @title Function to compute Bayes factors for general designs
##' @param formula a formula containing the full model for the analysis 
##'   (see Examples)
##' @param data a data frame containing data for all factors in the formula
##' @param whichRandom a character vector specifying which factors are random
##' @param whichModels which set of models to compare; see Details
##' @param neverExclude a character vector containing a regular expression (see
##' help for \link{regex} for details) that indicates which terms to always keep 
##' in the analysis
##' @param iterations How many Monte Carlo simulations to generate, if relevant
##' @param progress if \code{TRUE}, show progress with a text progress bar
##' @param rscaleFixed prior scale for standardized, reduced fixed effects. A 
##'   number of preset values can be given as strings; see Details.
##' @param rscaleRandom prior scale for standardized random effects
##' @param rscaleCont prior scale for standardized slopes
##' @param multicore if \code{TRUE} use multiple cores through the \code{doMC} 
##'   package. Unavailable on Windows.
##' @param method approximation method, if needed. See \code{\link{nWayAOV}} for
##'   details.
##' @param noSample if \code{TRUE}, do not sample, instead returning NA.
##' @return An object of class \code{BFBayesFactor}, containing the computed 
##'   model comparisons
##' @param callback callback function for third-party interfaces 
##' @author Richard D. Morey (\email{richarddmorey@@gmail.com})
##' @export
##' @references
##'   Rouder, J. N., Morey, R. D., Speckman, P. L., Province, J. M., (2012) 
##'   Default Bayes Factors for ANOVA Designs. Journal of Mathematical 
##'   Psychology.  56.  p. 356-374.
##'   
##'   Liang, F. and Paulo, R. and Molina, G. and Clyde, M. A. and 
##'   Berger, J. O. (2008). Mixtures of g-priors for Bayesian Variable 
##'   Selection. Journal of the American Statistical Association, 103, pp. 
##'   410-423
##'   
##' @note The function \code{generalTestBF} can compute Bayes factors for all 
##'   restrictions of a full model against the null 
##'   hypothesis that all effects are 0. The total number of tests 
##'   computed -- if all tests are requested -- will be \eqn{2^K-1}{2^K - 1} 
##'   for \eqn{K} factors or covariates.
##'   This number increases very quickly with the number of tested predictors. An option is included to
##'   prevent testing too many models: \code{options('BFMaxModels')}, which defaults to 50,000, is
##'   the maximum number of models that will be analyzed at once. This can
##'   be increased by increased using \code{\link{options}}.
##'   
##'   It is possible to reduce the number of models tested by only testing the 
##'   most complex model and every restriction that can be formed by removing 
##'   one factor or interaction using the \code{whichModels} argument. See the 
##'   help for \code{\link{anovaBF}} for details.
##'   
##' @examples
##' ## Puzzles example: see ?puzzles and ?anovaBF
##' data(puzzles)
##' ## neverExclude argument makes sure that participant factor ID
##' ## is in all models
##' result = generalTestBF(RT ~ shape*color + ID, data = puzzles, whichRandom = "ID", 
##' neverExclude="ID", progress=FALSE)
##' result
##' 
##' @keywords htest
##' @seealso \code{\link{lmBF}}, for testing specific models, and 
##'   \code{\link{regressionBF}} and \code{anovaBF} for other functions for
##'   testing multiple models simultaneously.
generalTestBF <- 
  function(formula, data, whichRandom = NULL, 
           whichModels = "withmain", neverExclude=NULL, iterations = 10000, progress = options()$BFprogress,
           rscaleFixed = "medium", rscaleRandom = "nuisance", rscaleCont="medium", multicore = FALSE, method="auto",noSample=FALSE, callback=function(...) as.integer(0))
  {
    checkFormula(formula, data, analysis = "lm")
    # pare whichRandom down to terms that appear in the formula
    whichRandom <- whichRandom[whichRandom %in% fmlaFactors(formula, data)[-1]]
    
    dataTypes <- createDataTypes(formula, whichRandom, data, analysis = "lm")

    models = enumerateGeneralModels(formula, whichModels, neverExclude, 
                                    includeBottom = whichModels!="top")
    
    if(length(models)>options()$BFMaxModels) stop("Maximum number of models exceeded (", 
                                                  length(models), " > ",options()$BFMaxModels ,"). ",
                                                  "The maximum can be increased by changing ",
                                                  "options('BFMaxModels').")
    
    if(multicore){
      message("Note: Progress bars and callbacks are suppressed when running multicore.")
      if(!requireNamespace("doMC", quietly = TRUE)){
        stop("Required package (doMC) missing for multicore functionality.")
      } 
      
      doMC::registerDoMC()
      if(foreach::getDoParWorkers()==1){
        warning("Multicore specified, but only using 1 core. Set options(cores) to something >1.")
      }
      
      bfs <- foreach::"%dopar%"(
        foreach::foreach(gIndex=models, .options.multicore=mcoptions),
        lmBF(gIndex,data = data, whichRandom = whichRandom, 
             rscaleFixed = rscaleFixed, rscaleRandom = rscaleRandom,
             rscaleCont = rscaleCont, iterations = iterations, method=method,
             progress=FALSE,noSample=noSample)
        )
    }else{ # Single core
      checkCallback(callback,as.integer(0))
      bfs = NULL
      myCallback <- function(prgs){
        frac <- (i - 1 + prgs/1000)/length(models)
        ret <- callback(frac*1000)
        return(as.integer(ret))
      }
      if(progress){
        pb = txtProgressBar(min = 0, max = length(models), style = 3)
      }else{
        pb = NULL
      }
      for(i in 1:length(models)){
        oneModel <- lmBF(models[[i]],data = data, whichRandom = whichRandom,
          rscaleFixed = rscaleFixed, rscaleRandom = rscaleRandom,
          rscaleCont = rscaleCont, iterations = iterations, 
          progress = FALSE, method = method,noSample=noSample,callback=myCallback)
        if(inherits(pb,"txtProgressBar")) setTxtProgressBar(pb, i)
        bfs = c(bfs,oneModel)
      }
      if(inherits(pb,"txtProgressBar")) close(pb)
      checkCallback(callback,as.integer(1000))
    }
    
    # combine all the Bayes factors into one BFBayesFactor object
    bfObj = do.call("c", bfs)
    
    if(whichModels=="top") bfObj = BFBayesFactorTop(bfObj)
    
    return(bfObj)
  }
