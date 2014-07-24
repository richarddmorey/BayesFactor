##' This function computes Bayes factors for contingency tables.
##' 
##' The Bayes factor provided by \code{contingencyTableBF} 
##' 
##' @title Function for Bayesian analysis of one- and two-sample designs
##' @param x an m by n matrix of counts (integers m,n > 1)
##' @param sampleType the sampling plan (see details)
##' @param priorConcentration prior concentration parameter (see details)
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
##' @author Tahira Jamil (\email{tahjamil@@gmail.com})
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
##' @note This is a note.
##' @examples
##' ## Fathers and Sons example
##' @seealso \code{\link{integrate}}, \code{\link{t.test}}

contingencyTableBF(x, sampleType, priorConcentration = 1, posterior = FALSE, ...)
{
  
  x = as.matrix(as.integer(x))

  numerator = switch(sampleType,
         poisson = BFcontingencyTable(type = "contingency table, poisson", 
                                               identifier = list(formula = "non-independence"), 
                                               prior=list(a=priorConcentration),
                                               shortName = "Non-indep.",
                                               longName = "Alternative, non-independence"),
         jointMulti = BFcontingencyTable(type = "contingency table, joint multinomial", 
                                         identifier = list(formula = "non-independence"), 
                                         prior=list(a=priorConcentration),
                                         shortName = "Non-indep.",
                                         longName = "Alternative, non-independence"),
         indepMulti = BFcontingencyTable(type = "contingency table, independent multinomial", 
                                         identifier = list(formula = "non-independence"), 
                                         prior=list(a=priorConcentration),
                                         shortName = "Non-indep.",
                                         longName = "Alternative, non-independence"),
         hypergeom = BFcontingencyTable(type = "contingency table, hypergeometric", 
                                        identifier = list(formula = "non-independence"), 
                                        prior=list(a=priorConcentration),
                                        shortName = "Non-indep.",
                                        longName = "Alternative, non-independence"),
         stop("Unknown value of sampleType (see help for contingencyBF).")
    )

    if(posterior){
      chains = posterior(numerator, data = x, ...)
      return(chains)
    }else{
      bf = compare(numerator = numerator, data = x)
      return(bf)
    }
}

