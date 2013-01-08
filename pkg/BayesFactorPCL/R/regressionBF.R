
enumerateRegressionModels = function(fmla, whichModels){
  trms <- attr(terms(fmla), "term.labels")
  ntrms <- length(trms)
  dv = deparse(fmla[[2]])
  if(ntrms == 1 ) whichModels = "all"
  
  if(whichModels=="top"){
    lst = combn2( trms, ntrms - 1 )
  }else if(whichModels=="all"){
    lst = combn2( trms, 1 )
  }else{
    stop("Unknown whichModels value: ",whichModels)
  }
  strng <- lapply(lst,function(el){
    paste(el,collapse=" + ")
  })
  fmla <- lapply(strng, function(el){
    formula(paste(dv,"~", el))
  })
  return(fmla)
}

createFullRegressionModel <- function(formula){
  factors = fmlaFactors(formula)[-1]
 
  dv = deparse(formula[[2]])
  
  RHS = paste(factors,collapse=" + ")
  strng = paste(dv, " ~ ", RHS, collapse = "")
  return(formula(strng))
}


regressionBF <- function(formula, data, whichModels = "all", progress=TRUE, rscaleCont = 1)
{
  checkFormula(formula, data, analysis = "regression")
  dataTypes <- createDataTypes(formula, whichRandom=c(), data, analysis = "regression")
  fmla <- createFullRegressionModel(formula)
  
  models <- enumerateRegressionModels(fmla, whichModels)
  
  if(progress) lapply = pblapply
  bfs <- lapply(models, lmBF, data = data, dataTypes = dataTypes,
                rscaleCont = rscaleCont)
  
  do.call("c", bfs)
}
