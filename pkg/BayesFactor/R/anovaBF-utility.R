
enumerateAnovaModelsWithMain<-function(factors){
  nFactors = length(factors)
  mb = monotoneBooleanNice(nFactors)
  mb = mb[-c(1,2),-ncol(mb)]
  myModels = apply(mb,1,function(v) rev(((length(v):1))[v]))
  myTerms = sapply(1:(2^nFactors-1),makeTerm,factors=factors)
  lapply(myModels, function(v) myTerms[v])
}

enumerateAnovaModels = function(fmla, whichModels, data){
  trms <- attr(terms(fmla, data = data), "term.labels")
  ntrms <- length(trms)
  dv = stringFromFormula(fmla[[2]])
  if(ntrms == 1 ) whichModels = "all"

  if(whichModels=="top"){
    lst = combn2( trms, ntrms - 1 )
  }else if(whichModels=='bottom'){
    lst = as.list(combn( trms, 1 ))
  }else if(whichModels=="all"){
    lst = combn2( trms, 1 )
  }else if(whichModels=="withmain"){
    lst = enumerateAnovaModelsWithMain( fmlaFactors(fmla, data)[-1] )
  }else{
    stop("Unknown whichModels value: ",whichModels)
  }
  strng <- sapply(lst,function(el){
    paste(el,collapse=" + ")
  })
  strng <- unique(strng)
  fmla <- lapply(strng, function(el){
    formula(paste(dv,"~", el))
  })
  return(fmla)
}

createFixedAnovaModel <- function(dataTypes, formula){
  fixedFactors <- names(dataTypes[dataTypes=="fixed"])
  fixedPart <- paste(fixedFactors,collapse="*")

  # get LHS of formula
  dv = stringFromFormula(formula[[2]])

  formula(paste(dv, "~", fixedPart, collapse=""))
}

addRandomModelPart <- function(formula, dataTypes, null = FALSE){
  randomFactors <- names(dataTypes[dataTypes=="random"])
  randomPart <- paste(randomFactors,collapse="+")

  fmla = stringFromFormula(formula)
  dv = stringFromFormula(formula[[2]])
  if(null){
    ret = formula(paste(dv, "~", randomPart, collapse=""))
  }else{
    ret = formula(paste(fmla, "+", randomPart, collapse=""))
  }
  return(ret)
}

