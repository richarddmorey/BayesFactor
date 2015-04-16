#' Compare two objects to see if they are the 'same', for some loose definition
#' of same
#' @param x first object
#' @param y second object
#' @return Returns \code{TRUE} or \code{FALSE}
setGeneric("%same%", function(x, y) standardGeneric("%same%"))

#' Find a model term in a vector of model terms
#' @param x the terms to be matched
#' @param table the terms to be matched against
#' @return A logical vector of the same length as x, indicating if a 
#' match was located for each element of x.
setGeneric("%termin%", function(x, table) standardGeneric("%termin%"))

#' Compare two models, with respect to some data
#' 
#' This method is used primarily in the backend, and will only rarely be called
#' by the end user. But see the examples below for a demonstration.
#' @param numerator first model
#' @param denominator second model (if omitted, compare to predefined null)
#' @param data data for the comparison
#' @param ... arguments passed to and from related methods
#' @return The compare function will return a model comparison object, typically
#'   a Bayes factor
#' @export
#' @docType methods
#' @rdname compare-methods
#' @aliases compare,BFoneSample,missing,data.frame-method 
#'  compare,BFlinearModel,BFlinearModel,data.frame-method
#'  compare,BFindepSample,missing,data.frame-method
#'  compare,BFlinearModel,missing,data.frame-method
#'  compare,BFmetat,missing,data.frame-method
#'  compare,BFproportion,missing,data.frame-method
#'  compare,BFcontingencyTable,BFcontingencyTable,data.frame-method
#'  compare,BFcontingencyTable,missing,data.frame-method
#'  compare,BFmcmc,BFmcmc,ANY-method
#'  compare,BFmcmc,missing,ANY-method
#' @examples
#' ## Sample from the posteriors for two models
#' data(puzzles)
#' 
#' ## Main effects model; result is a BFmcmc object, inheriting 
#' ## mcmc from the coda package
#' mod1 = lmBF(RT ~ shape + color + ID, data = puzzles, whichRandom = "ID", 
#'    progress = FALSE, posterior = TRUE, iterations = 1000)
#' 
#' plot(mod1)
#' 
#' ## Full model
#' mod2 = lmBF(RT ~ shape*color + ID, data = puzzles, whichRandom = "ID", 
#'    progress = FALSE, posterior = TRUE, iterations = 1000)
#' 
#' ## Each BFmcmc object contains the model used to generate it, so we
#' ## can compare them (data is not needed, it is contained in the objects):
#' 
#' compare(mod1, mod2) 
setGeneric("compare", function(numerator, denominator, data, ...) standardGeneric("compare"))

#' Recompute a Bayes factor computation or MCMC object. 
#' 
#' Take an object and redo the computation (useful for sampling). In cases where sampling is 
#' used to compute the Bayes factor, the estimate of the precision of new samples will be added
#' to the estimate precision of the old sample will be added to produce a new estimate of the 
#' precision.
#' @param x object to recompute
#' @param progress report progress of the computation?
#' @param multicore Use multicore, if available
#' @param callback callback function for third-party interfaces 
#' @param ... arguments passed to and from related methods
#' @return Returns an object of the same type, after repeating the sampling (perhaps with more iterations)
#' @export
#' @docType methods
#' @rdname recompute-methods
#' @examples
#' ## Sample from the posteriors for two models
#' data(puzzles)
#' 
#' ## Main effects model; result is a BFmcmc object, inheriting 
#' ## mcmc from the coda package
#' bf = lmBF(RT ~ shape + color + ID, data = puzzles, whichRandom = "ID", 
#'    progress = FALSE)
#'  
#' ## recompute Bayes factor object  
#' recompute(bf, iterations = 1000, progress = FALSE)
#'    
#' ## Sample from posterior distribution of model above, and recompute:
#' chains = posterior(bf, iterations = 1000, progress = FALSE)
#' newChains = recompute(chains, iterations = 1000, progress=FALSE)     
setGeneric("recompute", function(x, progress=options()$BFprogress, multicore = FALSE, callback = function(...) as.integer(0), ...) standardGeneric("recompute"))

#' Sample from the posterior distribution of one of several models.
#' 
#' This function samples from the posterior distribution of a \code{BFmodel}, 
#' which can be obtained from a \code{BFBayesFactor} object. If there is more 
#' than one numerator in the \code{BFBayesFactor} object, the \code{index} 
#' argument can be passed to select one numerator.
#' 
#' The data argument is used internally, and will y not be needed by 
#' end-users.
#' 
#' Note that if there are fixed effects in the model, the reduced 
#' parameterzation used internally (see help for \code{\link{anovaBF}}) is 
#' unreduced. For a factor with two levels, the chain will contain two effect 
#' estimates that sum to 0.
#' 
#' Two useful arguments that can be passed to related methods are \code{thin} 
#' and \code{columnFilter}, currently implemented for methods using 
#' \code{nWayAOV} (models with more than one categorical covariate, or a mix of 
#' categorical and continuous covariates). \code{thin}, an integer, will keep 
#' only every \code{thin} iterations. The default is \code{thin=1}, which keeps 
#' all iterations. Argument \code{columnFilter} is either \code{NULL} (for no 
#' filtering) or a character vector of extended regular expressions (see
#' \link{regex} help for details). Any column from an effect that matches one of
#' the filters will not be saved.
#' @param model or set of models from which to sample
#' @param index the index within the set of models giving the desired model
#' @param data the data to be conditioned on
#' @param iterations the number of iterations to sample
#' @param ... arguments passed to and from related methods
#' @return Returns an object containing samples from the posterior distribution 
#'   of the specified model
#' @export
#' @docType methods
#' @rdname posterior-methods
#' @examples
#' ## Sample from the posteriors for two models
#' data(sleep)
#' 
#' bf = lmBF(extra ~ group + ID, data = sleep, whichRandom="ID", progress=FALSE)  
#'    
#' ## sample from the posterior of the numerator model
#' ## data argument not needed - it is included in the Bayes factor object 
#' chains = posterior(bf, iterations = 1000, progress = FALSE)
#' 
#' plot(chains)
#' 
#' ## demonstrate column filtering by filtering out participant effects
#' data(puzzles)
#' bf = lmBF(RT ~ shape + color + shape:color + ID, data=puzzles)
#' chains = posterior(bf, iterations = 1000, progress = FALSE, columnFilter="^ID$")
#' colnames(chains) # Contains no participant effects
setGeneric("posterior", function(model, index, data, iterations, ...) standardGeneric("posterior"))

#' Extract the Bayes factor from an object
#' @param x object from which to extract the Bayes factors
#' @param logbf return the logarithm of the Bayes factors
#' @param onlybf return a vector of only the Bayes factors
#' @return Returns an object containing Bayes factors extracted from the object
#' @export
#' @docType methods
#' @rdname extractBF-methods
#' @examples
#' ## Sample from the posteriors for two models
#' data(puzzles)
#' 
#' bf = lmBF(RT ~ shape*color + ID, data = puzzles, whichRandom="ID", progress=FALSE)  
#'    
#' extractBF(bf)  
setGeneric("extractBF", function(x, logbf=FALSE, onlybf=FALSE) standardGeneric("extractBF"))

#' Extract the odds from an object
#' @param x object from which to extract
#' @param logodds return the logarithm
#' @param onlyodds return a vector of only the odds
#' @return Returns an object containing odds extracted from the object
#' @export
#' @docType methods
#' @rdname extractOdds-methods
setGeneric("extractOdds", function(x, logodds=FALSE, onlyodds=FALSE) standardGeneric("extractOdds"))

#' Extract the probabilities from an object
#' @param x object from which to extract
#' @param logprobs return the logarithm
#' @param onlyprobs return a vector of only the probabilities
#' @return Returns an object containing probabilities extracted from the object
#' @export
#' @docType methods
#' @rdname extractProbabilities-methods
setGeneric("extractProbabilities", function(x, logprobs=FALSE, onlyprobs=FALSE) standardGeneric("extractProbabilities"))

#' Set prior odds in an object
#' @docType methods
#' @rdname priorOdds-method
#' @param object object in which to set odds
#' @param value odds
setGeneric("priorOdds<-", function(object, value) standardGeneric("priorOdds<-"))

#' Set prior log odds in an object
#' @docType methods
#' @rdname priorLogodds-method
#' @param object object in which to set log odds
#' @param value log odds
setGeneric("priorLogodds<-", function(object, value) standardGeneric("priorLogodds<-"))

#' Filter the elements of an object according to some pre-specified criteria
#' @param x object
#' @param name regular expression to search name
#' @param perl logical. Should perl-compatible regexps be used? See ?grepl for details.
#' @param fixed logical. If TRUE, pattern is a string to be matched as is. See ?grepl for details.
#' @param ... arguments passed to and from related methods
#' @return Returns a filtered object
setGeneric("filterBF", function(x, name, perl = FALSE, fixed = FALSE, ...) standardGeneric("filterBF"))

