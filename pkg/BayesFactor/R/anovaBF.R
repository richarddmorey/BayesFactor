
##' This function computes Bayes factors for all main-effects and interaction 
##' contrasts in an ANOVA design.
##' 
##' Models, priors, and methods of computation are provided in Rouder et al. 
##' (2012).
##' 
##' The ANOVA model for a vector of observations \eqn{y} is \deqn{ y = \mu + X_1
##' \theta_1 + \ldots + X_p\theta_p +\epsilon,} where 
##' \eqn{\theta_1,\ldots,\theta_p} are vectors of main-effect and interaction 
##' effects, \eqn{X_1,\ldots,X_p} are corresponding design matrices, and 
##' \eqn{\epsilon} is a vector of zero-centered noise terms with variance 
##' \eqn{\sigma^2}.  Zellner and Siow (1980) inspired g-priors are placed on 
##' effects, but with a separate g-prior parameter for each covariate: 
##' \deqn{\theta_1~N(0,g_1\sigma^2), \ldots,  \theta_p~N(0,g_p \sigma^2).}  A 
##' Jeffries prior is placed on \eqn{\mu} and \eqn{\sigma^2}.  Independent 
##' scaled inverse-chi-square priors with one degree of freedom are placed on 
##' \eqn{g_1,\ldots,g_p}.  The square-root of the scale for g's corresponding to
##' fixed and random effects is given by \code{rscaleFixed} and 
##' \code{rscaleRandom}, respectively.
##' 
##' When a factor is treated as random, there are as many main effect terms in 
##' the vector \eqn{\theta} as levels.  When a factor is treated as fixed, the 
##' sums-to-zero linear constraint is enforced by centering the corresponding 
##' design matrix, and there is one fewer main effect terms as levels.  The 
##' Cornfield-Tukey model of interactions is assumed.  Details are provided in 
##' Rouder et al. (2012)
##' 
##' Bayes factors are computed by integrating the likelihood with respect to the
##' priors on parameters.  The integration of all parameters except 
##' \eqn{g_1,\ldots,g_p} may be expressed in closed-form; the integration of 
##' \eqn{g_1,\ldots,g_p} is performed through Monte Carlo sampling, and 
##' \code{iterations} is the number of iterations used to estimate the Bayes 
##' factor.
##' 
##' \code{anovaBF} computes Bayes factors for either all submodels or select 
##' submodels missing a single main effect or covariate, depending on the 
##' argument \code{whichModels}. If no random factors are specified, the null 
##' model assumed by \code{anovaBF} is the grand-mean only model. If random 
##' factors are specified, the null model is the model with an additive model on
##' all random factors, plus a grand mean. Thus, \code{anovaBF} does not 
##' currently test random factors. Testing random factors is possible with 
##' \code{\link{lmBF}}.
##' 
##' The argument \code{whichModels} controls which models are tested. Possible 
##' values are 'all', 'withmain', 'top', and 'bottom'. Setting 
##' \code{whichModels} to 'all' will test all models that can be created by 
##' including or not including a main effect or interaction. 'top' will test all
##' models that can be created by removing or leaving in a main effect or 
##' interaction term from the full model. 'bottom' creates models by adding 
##' single factors or interactions to the null model. 'withmain' will test all 
##' models, with the constraint that if an interaction is included, the 
##' corresponding main effects are also included.
##' 
##' For the \code{rscaleFixed} and \code{rscaleRandom} arguments, several named
##' values are recognized: "medium", "wide", and "ultrawide", corresponding to
##' \eqn{r} scale values of 1/2, \eqn{\sqrt{2}/2}{sqrt(2)/2}, and 1,
##' respectively. In addition, \code{rscaleRandom} can be set to the "nuisance",
##' which sets \eqn{r=1} (and is thus equivalent to "ultrawide"). The "nuisance"
##' setting is for medium-to-large-sized effects assumed to be in the data but 
##' typically not of interest, such as variance due to participants.
##' @title Function to compute Bayes factors for ANOVA designs
##' @param formula a formula containing all factors to include in the analysis 
##'   (see Examples)
##' @param data a data frame containing data for all factors in the formula
##' @param whichRandom a character vector specifying which factors are random
##' @param whichModels which set of models to compare; see Details
##' @param iterations How many Monte Carlo simulations to generate, if relevant
##' @param progress if \code{TRUE}, show progress with a text progress bar
##' @param rscaleFixed prior scale for standardized, reduced fixed effects. A 
##'   number of preset values can be given as strings; see Details.
##' @param rscaleRandom prior scale for standardized random effects
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
##' @references Gelman, A. (2005) Analysis of Variance---why it is more 
##'   important than ever.  Annals of Statistics, 33, pp. 1-53.
##'   
##'   Rouder, J. N., Morey, R. D., Speckman, P. L., Province, J. M., (2012) 
##'   Default Bayes Factors for ANOVA Designs. Journal of Mathematical 
##'   Psychology.  56.  p. 356-374.
##'   
##'   Zellner, A. and Siow, A., (1980) Posterior Odds Ratios for Selected 
##'   Regression Hypotheses.  In Bayesian Statistics: Proceedings of the First 
##'   Interanational Meeting held in Valencia (Spain).  Bernardo, J. M., 
##'   Lindley, D. V., and Smith A. F. M. (eds), pp. 585-603.  University of 
##'   Valencia.
##'   
##' @note The function \code{anovaBF} will compute Bayes factors for all 
##'   possible combinations of fixed factors and interactions, against the null 
##'   hypothesis that \emph{all} effects are 0. The total number of tests 
##'   computed will be \eqn{2^{2^K - 1}}{2^(2^K - 1)} for \eqn{K} fixed factors.
##'   This number increases very quickly with the number of factors. For 
##'   instance, for a five-way ANOVA, the total number of tests exceeds two 
##'   billion. Even though each test takes a fraction of a second, the time 
##'   taken for all tests could exceed your lifetime. An option is included to
##'   prevent this: \code{options('BFMaxModels')}, which defaults to 50,000, is
##'   the maximum number of models that `anovaBF` will analyze at once. This can
##'   be increased by increasing the option value.
##'   
##'   It is possible to reduce the number of models tested by only testing the 
##'   most complex model and every restriction that can be formed by removing 
##'   one factor or interaction using the \code{whichModels} argument. Setting 
##'   this argument to 'top' reduces the number of tests to \eqn{2^K-1}, which 
##'   is more manageable. The Bayes factor for each restriction against the most
##'   complex model can be interpreted as a test of the removed 
##'   factor/interaction. Setting \code{whichModels} to 'withmain' will not 
##'   reduce the number of tests as much as 'top' but the results may be more 
##'   interpretable, since an interaction is only allowed when all interacting 
##'   effects (main or interaction) are also included in the model.
##'   
##' @examples
##' ## Classical example, taken from t.test() example
##' ## Student's sleep data
##' data(sleep)
##' plot(extra ~ group, data = sleep)
##' 
##' ## traditional ANOVA gives a p value of 0.00283
##' summary(aov(extra ~ group + Error(ID/group), data = sleep))
##' 
##' ## Gives a Bayes factor of about 11.6
##' ## in favor of the alternative hypothesis
##' anovaBF(extra ~ group + ID, data = sleep, whichRandom = "ID", 
##'     progress=FALSE)
##' 
##' ## Demonstrate top-down testing
##' data(puzzles)
##' result = anovaBF(RT ~ shape*color + ID, data = puzzles, whichRandom = "ID", 
##'     whichModels = 'top', progress=FALSE)
##' result
##' 
##' ## In orthogonal designs, the top down Bayes factor can be
##' ## interpreted as a test of the omitted effect
##' @keywords htest
##' @seealso \code{\link{lmBF}}, for testing specific models, and 
##'   \code{\link{regressionBF}} for the function similar to \code{anovaBF} for 
##'   linear regression models.
anovaBF <- 
  function(formula, data, whichRandom = NULL, 
           whichModels = "withmain", iterations = 10000, progress = options()$BFprogress,
           rscaleFixed = "medium", rscaleRandom = "nuisance", multicore = FALSE, method="auto", noSample=FALSE, callback=function(...) as.integer(0))
  {
    checkFormula(formula, data, analysis = "anova")
    # pare whichRandom down to terms that appear in the formula
    whichRandom <- whichRandom[whichRandom %in% fmlaFactors(formula, data)[-1]]
    if(all(fmlaFactors(formula, data)[-1] %in% whichRandom)){
      # No fixed factors!
      bf = lmBF(formula, data, whichRandom,rscaleFixed,rscaleRandom, progress=progress, method=method,noSample=noSample)  
      return(bf)
    }
    
    dataTypes <- createDataTypes(formula, whichRandom, data, analysis = "anova")
    fmla <- createFixedAnovaModel(dataTypes, formula)
    
    models <- enumerateAnovaModels(fmla, whichModels, data)
    
    if(length(models)>options()$BFMaxModels) stop("Maximum number of models exceeded (", 
                                                  length(models), " > ",options()$BFMaxModels ,"). ",
                                                  "The maximum can be increased by changing ",
                                                  "options('BFMaxModels').")
    
    if(length(whichRandom) > 0 ){
      models <- lapply(models, addRandomModelPart, dataTypes = dataTypes)
      models <- c(models, addRandomModelPart(fmla, dataTypes, null=TRUE))
    }
    
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
             iterations = iterations, method=method,progress=FALSE,noSample=noSample)
        )
      
    }else{ # Single core

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
                         iterations = iterations, progress = FALSE, method = method,
                         noSample=noSample,callback=myCallback)
        if(inherits(pb,"txtProgressBar")) setTxtProgressBar(pb, i)
        bfs = c(bfs,oneModel)
      }
      if(inherits(pb,"txtProgressBar")) close(pb)      
      
    }
    
    # combine all the Bayes factors into one BFBayesFactor object
    bfObj = do.call("c", bfs)
    # If we have random effects, make those the denominator
    if(length(whichRandom) > 0) bfObj = bfObj[-length(bfObj)] / bfObj[length(bfObj)]
    
    if(whichModels=="top") bfObj = BFBayesFactorTop(bfObj)
    
    return(bfObj)
  }


