##' This function computes Bayes factors for contingency tables.
##' 
##' The Bayes factor provided by \code{contingencyTableBF} tests the independence assumption in 
##' contingency tables under various sampling plans.
##' 
##' @title Function for Bayesian analysis of one- and two-sample designs
##' @param x an m by n matrix of counts (integers m,n > 1)
##' @param sampleType the sampling plan (see details)
##' @param fixedMargin for the independent multinomial sampling plan, which margin is fixed ("rows" or "cols")
##' @param priorConcentration prior concentration parameter (see details)
##' @param posterior if \code{TRUE}, return samples from the posterior instead 
##'   of Bayes factor
##' @param callback callback function for third-party interfaces 
##' @param ... further arguments to be passed to or from methods.
##' @return If \code{posterior} is \code{FALSE}, an object of class 
##'   \code{BFBayesFactor} containing the computed model comparisons is 
##'   returned. 
##'   
##'   If \code{posterior} is \code{TRUE}, an object of class \code{BFmcmc},
##'   containing MCMC samples from the posterior is returned.
##' @export
##' @keywords internal htest
##' @author Richard D. Morey (\email{richarddmorey@@gmail.com})
##' @author Tahira Jamil (\email{tahjamil@@gmail.com})
##' @references Gunel, E. and Dickey, J., (1974) 
##' Bayes Factors for Independence in Contingency Tables. Biometrika, 61, 545-557
##' 
##' @note This will be implemented in a future version.
##' @examples
##' \dontrun{
##' 
##' data<-matrix(c(10,3,2,15),c(2,2))
##' 
##' ## Assume poisson sampling scheme
##' contingencyTableBF(data, "poisson")
##' }

contingencyTableBF <- function(x, sampleType, fixedMargin = NULL, priorConcentration = 1, posterior = FALSE,  callback = function(...) as.integer(0), ...)
{
  
  x.mat = as.matrix(as.integer(x))
  dim(x.mat) = dim(x)
  x = as.data.frame(x.mat)
  
  if( sampleType == "indepMulti" )
    if( is.null(fixedMargin) ){
      stop("Argument fixedMargin ('rows' or 'cols') required with independent multinomial sampling plan.")
    }else if( !(fixedMargin %in% c("rows","cols")) ){
      stop("Argument fixedMargin must be either 'rows' or 'cols'.")
    }

  
  numerator = switch(sampleType,
         poisson = BFcontingencyTable(type = "poisson", 
                                               identifier = list(formula = "non-independence"), 
                                               prior=list(a=priorConcentration),
                                               shortName = paste0("Non-indep. (a=",priorConcentration,")"),
                                               longName = paste0("Alternative, non-independence, a = ", priorConcentration)),
         jointMulti = BFcontingencyTable(type = "joint multinomial", 
                                         identifier = list(formula = "non-independence"), 
                                         prior=list(a=priorConcentration),
                                         shortName = paste0("Non-indep. (a=",priorConcentration,")"),
                                         longName = paste0("Alternative, non-independence, a = ", priorConcentration)),
         indepMulti = BFcontingencyTable(type = "independent multinomial", 
                                         identifier = list(formula = "non-independence", fixedMargin = fixedMargin), 
                                         prior=list(a=priorConcentration, fixedMargin = fixedMargin),
                                         shortName = paste0("Non-indep. (a=",priorConcentration,")"),
                                         longName = paste0("Alternative, non-independence, a = ", priorConcentration)),
         hypergeom = BFcontingencyTable(type = "hypergeometric", 
                                        identifier = list(formula = "non-independence"), 
                                        prior=list(a=priorConcentration),
                                        shortName = paste0("Non-indep. (a=",priorConcentration,")"),
                                        longName = paste0("Alternative, non-independence, a = ", priorConcentration)),
         stop("Unknown value of sampleType (see help for contingencyTableBF).")
    )

    if(posterior){
      chains = posterior(numerator, data = x, callback = callback, ...)
      return(chains)
    }else{
      bf = compare(numerator = numerator, data = x)
      return(bf)
    }
}

