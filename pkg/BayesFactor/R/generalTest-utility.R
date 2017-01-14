requiredFor = Vectorize(function(t1,t2){
  t1 = unique(unlist(strsplit(t1,":",fixed=TRUE)))  
  t2 = unique(unlist(strsplit(t2,":",fixed=TRUE)))  
  return(all(t1 %in% t2))
},c("t1","t2"))

##' Generate lists of nested models, given a model formula
##' 
##' This is a backend function not intended for users. It is exposed for third-party
##' applications.
##' @title Function for generation of nested linear models
##' @param fmla formula for the "full" model
##' @param whichModels which subsets of models to generate 
##' @param neverExclude a character vector of terms to never remove 
##' @param includeBottom Include the base model containing only \code{neverExclude} terms
##' @param data a data frame containing the columns mentioned in \code{fmla}
##' @keywords internal
enumerateGeneralModels = function(fmla, whichModels, neverExclude=NULL, includeBottom=TRUE, data=NULL){
  trms <- attr(terms(fmla, data = data), "term.labels")
  
  # Remove everything we never exclude, to replace them later
  logicalToInclude = filterVectorLogical(neverExclude,trms)
  if(any(logicalToInclude)){
    alwaysIncluded = trms[logicalToInclude]
    if(whichModels=="withmain"){
      rq = matrix(outer(trms,alwaysIncluded,requiredFor),nrow=length(trms))
      rq = apply(rq,1,any)
      logicalToInclude[rq] = TRUE 
      alwaysIncluded = unique(c(trms[rq], alwaysIncluded))
    }
    alwaysIncludedString = paste(alwaysIncluded,collapse=" + ")
  }
  trms =  trms[!logicalToInclude]
  
  ntrms <- length(trms)
  dv = stringFromFormula(fmla[[2]])
  if(ntrms == 0 )  return(list(fmla))
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
  return(fmlaOut)
}


possibleRestrictionsWithMainGeneralFallback <- function(trms, alwaysKept){
  ntrms = length( trms )
  if(ntrms==1) return(NULL)
  thisLevelRestrictions = lapply(1:ntrms,function(i, trms, alwaysKept ){
    removed = unlist(strsplit(trms[i],":",fixed=TRUE))
    remaining = trms[-i]
    containsRemoved = sapply(remaining,function(el){
      splt = unlist(strsplit(el,":",fixed=TRUE))
      all(removed %in% splt)
    })    
    if(any(containsRemoved)){
      return(NULL)
    }else{
      return(remaining)
    }
  },trms=trms,alwaysKept=alwaysKept)
  
  thisLevelRestrictions = thisLevelRestrictions[!sapply(thisLevelRestrictions, is.null)]
  nextLevelRestrictions <- lapply(thisLevelRestrictions, possibleRestrictionsWithMainGeneralFallback, alwaysKept=alwaysKept)
  
  bothLevelRestrictions = c(thisLevelRestrictions,unlist(nextLevelRestrictions,recursive=FALSE))
  bothLevelRestrictions = bothLevelRestrictions[!sapply(bothLevelRestrictions, is.null)]
  return(unique(bothLevelRestrictions))                                                                  
}

possibleRestrictionsWithMainGeneral <- function(trms, alwaysKept=NULL){
  if(length(trms)==1) return(list(trms))
  myFactors = unique(unlist(strsplit(trms,":",fixed=TRUE)))  
  nFactors = length(myFactors)
  
  # If there are more than 5 factors involved then the total number of models is 
  # over 7 million; fall back to search-based method (becase there may not be
  # that many)
  if(nFactors>getOption('BFfactorsMax', 5)){
    warning("Falling back to slow recursive method of enumerating models due to many factors.")
    retList = possibleRestrictionsWithMainGeneralFallback(trms, alwaysKept)
    return(c(retList,list(trms)))
  }
  
  myTerms = sapply(1:(2^nFactors-1),makeTerm,factors=myFactors)
  
  # These terms MUST be in the model
  toKeep = rev(myTerms) %termin% alwaysKept
  
  # These terms should not be in the model, because they were not in the original
  # specification
  toDiscard = !(rev(myTerms) %termin% c(trms, alwaysKept))
  
  # The specified full model, specified as TRUE/FALSE on the list of terms
  row = rep(FALSE,2^nFactors-1)
  row[termMatch(c(trms,alwaysKept),rev(myTerms))]=TRUE
  
  # Get all possible models with nFactors factors, as matrix
  mb = monotoneBooleanNice(nFactors)
  mb = mb[-c(1,2),-ncol(mb)]
  
  if(dim(mb)[1]==0) stop("No models left in analysis. Please check that your model is valid under 'withmain'.")
  
  # Remove rows that have include
  invalidRows = apply(mb,1,function(v) any(v[toDiscard]))
  mb = matrix(mb[!invalidRows,], ncol = ncol(mb))
  
  if(dim(mb)[1]==0) stop("No models left in analysis. Please check that your model is valid under 'withmain'.")
  
  # Remove rows that do NOT include required terms
  validRows = apply(mb,1,function(v) all(v[toKeep]))
  mb = matrix(mb[validRows,], ncol = ncol(mb))

  if(dim(mb)[1]==0) stop("No models left in analysis. Please check that your model is valid under 'withmain'.")
  
  # Get all submodels of the specified model
  subMods = subModelsMatrix(row,mb)
  
  # Turn logicals into term numbers 
  myModels = apply(subMods,1,function(v) rev(((length(v):1))[v]))
  
  if(nrow(subMods)==1){
    retVec = myTerms[myModels]
    # eliminate anything that should be always kept, to be added in later
    retVec = retVec[!(retVec %termin% alwaysKept)]
    if(length(retVec)>0){
      return(list(retVec))
    }else{
      return(NULL)
    }
  }else{
    retList = lapply(myModels, function(v){
      retVec = myTerms[v]
      # eliminate anything that should be always kept, to be added in later
      retVec = retVec[!(retVec %termin% alwaysKept)]
      if(length(retVec)>0){
        return(retVec)
      }else{
        return(NULL)
      }
    })
    return(retList[!sapply(retList, is.null)])
  }
}

subModelsMatrix<-function(row,monoBool){
  if(length(row) != ncol(monoBool)) stop("Invalid number of terms in submodel")
  rows = apply(monoBool,1,function(v) all(row[v]))
  rowNum = which(apply(monoBool,1,function(v) all(row==v)))
  if(!any(rows)) return(matrix(row,nrow=1))
  retMatrix = matrix(monoBool[rows,],nrow=sum(rows))
  if(length(rowNum)>0){
    return(retMatrix)
  }else{
    return(rbind(retMatrix,row))
  }
}
