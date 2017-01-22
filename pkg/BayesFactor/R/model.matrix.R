
designMatrix = function(bf, ...){

  model = bf@numerator[[1]]

  if( class(model) != "BFlinearModel" ) stop("Model matrix not defined for this model type.")
  designMatrixJZS_LM(bf, ...)
}

designMatrixJZS_LM = function(bf, ...){

  model = bf@numerator[[1]]
  data = bf@data
  dataTypes = model@dataTypes

  formula = formula(model@identifier$formula)
  checkFormula(formula, data, analysis = "lm")

  factors = fmlaFactors(formula, data)[-1]
  nFactors = length(factors)

  if( nFactors == 0 ){
    X = matrix(1,nrow(data),1)
    gMap = c(intercept=NA)
  }else{
    # Remove "as.matrix" when sparse matrix support is added
    X = as.matrix(fullDesignMatrix(formula, data, dataTypes))
    gMap = createGMap(formula, data, dataTypes)
    X = cbind(1,X)
    gMap = c(intercept=NA, gMap)
  }

  attr(X,"gMap") = gMap
  return(X)
}

#' Design matrices for Bayes factor linear models analyses.
#'
#' This function returns the design matrix used for computation of the Bayes factor
#' for the numerator of a \code{BFBayesFactor} object. There must not be more
#' than one numerator in the \code{BFBayesFactor} object.
#' @param object a BayesFactor object with a single numerator
#' @param ... arguments passed to and from related methods
#' @return Returns the design matrix for the corresponding model. The 'gMap' attribute of the returned
#' matrix contains the mapping from columns of the design matrix to g parameters
#' @export
#' @docType methods
#' @rdname model.matrix-methods
#' @aliases model.matrix,BFBayesFactor
#' @references Rouder, J. N., Morey, R. D., Speckman, P. L., Province, J. M., (2012)
#'   Default Bayes Factors for ANOVA Designs. Journal of Mathematical
#'   Psychology.  56.  p. 356-374.
#' @examples
#' ## Gets the design matrix for a simple analysis
#' data(sleep)
#'
#' bf = anovaBF(extra ~ group + ID, data = sleep, whichRandom="ID", progress=FALSE)
#' X = model.matrix(bf)
#'
#' ## Show dimensions of X (should be 20 by 12)
#' dim(X)
setMethod('model.matrix', signature(object = "BFBayesFactor"),
          function(object, ...){
            if(length(object)>1) stop("Must specify single model.")
            designMatrix(object, ...)
          }
)

#' @rdname model.matrix-methods
#' @aliases model.matrix,BFBayesFactor
setMethod('model.matrix', signature(object = "BFBayesFactorTop"),
          function(object, ...){
            if(length(object)>1) stop("Must specify single model.")
            designMatrix(as.BFBayesFactor(object), ...)
          }
)
