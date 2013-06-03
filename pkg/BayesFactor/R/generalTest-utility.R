enumerateGeneralModels = function(fmla, whichModels, neverExclude=NULL,includeBottom=TRUE){
  trms <- attr(terms(fmla), "term.labels")
  
  # Remove everything we never exclude, to replace them later
  logicalToInclude = filterVectorLogical(neverExclude,trms)
  if(any(logicalToInclude)){
    alwaysIncluded = trms[logicalToInclude]
    alwaysIncludedString = paste(alwaysIncluded,collapse=" + ")
  }
  trms =  trms[!logicalToInclude]
  
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
    if(any(logicalToInclude)){
      lst = possibleRestrictionsWithMainGeneral( trms, alwaysIncluded )
    }else{
      lst = possibleRestrictionsWithMainGeneral( trms, NULL )
    }
  }else{
    stop("Unknown whichModels value: ",whichModels)
  }
  strng <- sapply(lst,function(el, suffix){
    paste(el,collapse=" + ")
  })
  # Add back in the terms to always include
  if(any(logicalToInclude)){
    strng <- sapply(strng,function(el, suffix){
      paste(el,suffix,collapse=" + ",sep=" + ")
    },suffix=alwaysIncludedString)
    if(includeBottom)
      strng <- c(strng,alwaysIncludedString)
  }  
  strng <- unique(strng)
  fmlaOut <- lapply(strng, function(el){
    formula(paste(dv,"~", el))
  })
  if(whichModels=="withmain"){
    return(c(fmla,fmlaOut))
  }else{
    return(fmlaOut)
  }
}


possibleRestrictionsWithMainGeneral <- function( trms, alwaysKept=NULL ){
  ntrms = length( trms )
  if(ntrms==1) return(NULL)
  thisLevelRestrictions = lapply(1:ntrms,function(i, trms, alwaysKept ){
    removed = trms[i]
    remaining = trms[-i]
    searchTerms = c(paste("^",removed,":",sep=""),
                    paste(":",removed,"$",sep=""),
                    paste(":",removed,":",sep="")
                    )
    containsRemoved = filterVectorLogical(searchTerms, c(remaining,alwaysKept))
    if(any(containsRemoved)){
      return(NULL)
    }else{
      return(remaining)
    }
  },trms=trms,alwaysKept=alwaysKept)
  
  thisLevelRestrictions = thisLevelRestrictions[!sapply(thisLevelRestrictions, is.null)]
  nextLevelRestrictions <- lapply(thisLevelRestrictions, possibleRestrictionsWithMainGeneral, alwaysKept=alwaysKept)

  bothLevelRestrictions = c(thisLevelRestrictions,unlist(nextLevelRestrictions,recursive=FALSE))
  bothLevelRestrictions = bothLevelRestrictions[!sapply(bothLevelRestrictions, is.null)]
  return(unique(bothLevelRestrictions))                                                                  
}

combine.formulas <- function(formula, with) {
  rh.formula <- paste(attr(terms(formula), "term.labels"), collapse = "+")
  rh.with <- paste(attr(terms(with), "term.labels"), collapse = "+")
  as.formula(paste(formula[[2]], "~", rh.formula, "+", rh.with))
}
