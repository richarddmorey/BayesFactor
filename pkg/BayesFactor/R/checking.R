createDataTypes <- function(formula, whichRandom, data, analysis){
  factors <- rownames(attr(terms(formula, data = data),"factors"))[-1]
  factors <- unlist(decomposeTerms(factors))
  cnames <- colnames(data)

  # check status of data columns
  types = sapply(cnames, function(name, data){
    ifelse(is.factor(data[,name]), "fixed", "continuous")
  }, data = data)

  # restrict to only columns of interest
  types = types[ names(types) %in% factors ]
  if(length(types) <1 ) return(c())

  if( any(types[ names(types) %in% whichRandom ] == "continuous") )
    stop("Nonfactors are specified as random.")
  if(length(whichRandom)>0)
    types[ names(types) %in% whichRandom ] = "random"

  #### check various analysis types

  ## ANOVA can only accept factors
  if( any(types=="continuous") & analysis == "anova" )
    stop("anovaBF() cannot be used with nonfactor independent variables. Use lmBF(), regressionBF(), or generalTestBF() instead.")

  ## regression can only accept nonfactors
  if( any(types %in% c("fixed", "random")) & analysis == "regression" )
    stop("regressionBF() cannot be used with factor independent variables. Use lmBF(), anovaBF(), or generalTestBF() instead.")

  #### End checking analysis types
  return(types)
}



checkFormula <- function(formula, data, analysis){
  if(length(formula) < 3) stop("LHS of formula must be given.")
  cnames = colnames(data)

  dv = stringFromFormula(formula[[2]])

  if(!is.numeric(data[,dv])) stop("Dependent variable must be numeric.")
  if(any(is.na(data[,dv])) | any(is.infinite(data[,dv]))) stop("Dependent variable must not contain missing or infinite values.")

  factors = fmlaFactors(formula, data)
  terms = colnames(attr(terms(formula, data = data),"factors"))
  decom = decomposeTerms(terms)
  terms = unlist(decom)

  vars = rownames(attr(terms(formula, data = data),"factors"))
  vars = unlist(decomposeTerms(vars))

  if(any(is.na(data[,vars]))) stop("Predictors must not contain missing values.")

  if(is.null(factors)) return()
  if(factors[1] %in% terms) stop("Dependent variable cannot be a predictor.")
  if(!all(factors %in% cnames)) stop("Some variables missing in data frame.")

  if(analysis=="regression"){
    RHS = stringFromFormula(formula[[3]])
    lengths = sapply(decom, length)
    if (any(lengths > 1)) stop("Interactions not allowed in regressionBF (try generalTestBF).")
  }

  if(analysis=="lm" | analysis=="anova" | analysis == "regression" | analysis == "indept")
    if(attr(terms(formula, data = data),"intercept") == 0) stop("Formula must include intercept.")

  if(analysis=="indept"){
    if( length(decom) > 1 ) stop("Indep. groups t test can only support 1 factor as predictor.")
    if(length(decom[[1]]) > 1) stop("Interaction terms are not allowed in t test.")
    if(nlevels(factor(data[,terms])) > 2) stop("Indep. groups t test requires a factor with exactly 2 levels.")
  }
  invisible()
}

checkEffects <- function(effects, data, dataTypes){
  if(!all(effects %in% colnames(data))) stop("Term in formula missing in data")
  if(!all(effects %in% names(dataTypes))) stop("Term in formula missing in dataTypes")
  # add more checking code here
  # most importantly, to check consistancy of data factors and dataTypes
  # no factors should be labeled as continuous, etc
}
