checkPossibleEffect <- function(effect,lowerLevel,struc)
{
  order = struc$nEffs[struc$sums==effect]
  if(order==1) return(TRUE)
  vec = struc$mat[struc$sums==effect,]
  vec = vec[vec != 0]
  necessaryEffects = combn(vec,order-1,sum)
  if(all(necessaryEffects %in% lowerLevel)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}  

enumeratePossibleModels <- function(lowerLevel=NULL,order,struc){
  nFac = struc$nFac
  if(order==1){
    mains = 2^(1:nFac - 1)
    combs = combn2(mains,1)
    nextLevel = lapply(combs,enumeratePossibleModels,order=2,struc=struc)
    nextLevel = nextLevel[!sapply(nextLevel, is.null)]
    nextLevel = unlist(nextLevel,recursive=FALSE)
    return(c(combs, nextLevel))
  }
  if(sum(struc$nEffs[match(lowerLevel,struc$sums)]==(order-1))<order) return()
  possibleEffects = struc$sums[struc$nEffs==order]
  checkEffects = sapply(possibleEffects,checkPossibleEffect,lowerLevel=lowerLevel,struc=struc)
  possibleEffects = possibleEffects[checkEffects]
  
  if(length(possibleEffects)==0) return()
  
  if(order==nFac){
    return(list(c(lowerLevel,possibleEffects)))
  }else{
    if(length(possibleEffects)>1){
      combs = combn2(possibleEffects,1)
    }else{
      combs = possibleEffects
    }
    names(combs)=NULL
    models = lapply(combs,function(effects,lowerLevel) c(lowerLevel,effects), lowerLevel=lowerLevel)
    nextLevel = lapply(models,enumeratePossibleModels,order=order+1,struc=struc)
    nextLevel = nextLevel[!sapply(nextLevel, is.null)]
    nextLevel = unlist(nextLevel,recursive=FALSE)
    return(c(models,nextLevel))
  }
}

# This still uses model numbers. I want to rework it, but I'm not sure how to easily
# do this with factor names. For now, I'll just convert them at the end
enumerateAnovaModelsWithMain<-function(factors){
  nFac = length(factors)
  effs = 2^(1:nFac - 1)
  struc = list()
  struc$mat = expand.grid(lapply(effs,function(el) c(0,el)))[-1,]
  struc$sums = rowSums(struc$mat)
  struc$nEffs = rowSums(struc$mat != 0)
  struc$nFac = nFac
  mods = enumeratePossibleModels(,1,struc)
  mods = lapply(mods, modelIndicatorToName, factors = factors)  
  return(mods)
}

modelIndicatorToName <- function(indic, factors){
  effects <- sapply(indic, 
                function(num, factors){
                  facs = factors[binary(num,dim=length(factors))$dicotomy]
                  paste(facs,collapse=":")
                }, factors = factors)
  return(effects)
}


enumerateAnovaModels = function(fmla, whichModels){
  trms <- attr(terms(fmla), "term.labels")
  ntrms <- length(trms)
  dv = deparse(fmla[[2]])
  if(ntrms == 1 ) whichModels = "all"
  
  if(whichModels=="top"){
    lst = combn2( trms, ntrms - 1 )
  }else if(whichModels=="all"){
    lst = combn2( trms, 1 )
  }else if(whichModels=="withmain"){
    lst = enumerateAnovaModelsWithMain( fmlaFactors(fmla)[-1] )
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
  dv = deparse(formula[[2]])
  
  formula(paste(dv, "~", fixedPart, collapse=""))
}

addRandomModelPart <- function(formula, dataTypes, null = FALSE){
  randomFactors <- names(dataTypes[dataTypes=="random"])
  randomPart <- paste(randomFactors,collapse="+")
  
  fmla = deparse(formula)
  dv = deparse(formula[[2]])
  if(null){
    ret = formula(paste(dv, "~", randomPart, collapse=""))
  }else{
    ret = formula(paste(fmla, "+", randomPart, collapse=""))    
  }
}


anovaBF <- 
  function(formula, data, whichRandom = NULL, 
           whichModels = "withmain", iterations = 10000, progress = TRUE,
           rscaleFixed = "medium", rscaleRandom = 1, multicore = FALSE)
  {
    checkFormula(formula, data, analysis = "anova")
    # pare whichRandom down to terms that appear in the formula
    whichRandom <- whichRandom[whichRandom %in% fmlaFactors(formula)[-1]]
    if(all(fmlaFactors(formula)[-1] %in% whichRandom)){
      # No fixed factors!
      bf = lmBF(formula, data, whichRandom,rscaleFixed,rscaleRandom)  
      return(bf)
    }
    
    dataTypes <- createDataTypes(formula, whichRandom, data, analysis = "anova")
    fmla <- createFixedAnovaModel(dataTypes, formula)
    
    models <- enumerateAnovaModels(fmla, whichModels)
    
    if(length(whichRandom) > 0 ){
      models <- lapply(models, addRandomModelPart, dataTypes = dataTypes)
      models <- c(models, addRandomModelPart(fmla, dataTypes, null=TRUE))
    }
    
    if(multicore){
      if(progress) warning("Progress bars are suppressed when running multicore.")
      if(!require(doMC)){
        stop("Required package (doMC) missing for multicore functionality.")
      } 
      
      registerDoMC()
      if(getDoParWorkers()==1){
        warning("Multicore specified, but only using 1 core. Set options(cores) to something >1.")
      }
      
      bfs <- foreach(gIndex=models, .options.multicore=mcoptions) %dopar% 
        lmBF(gIndex,data = data, whichRandom = whichRandom, 
             rscaleFixed = rscaleFixed, rscaleRandom = rscaleRandom,
             iterations = iterations)
      
    }else{ # Single core
      
      if(progress) lapply = pblapply
      bfs <- lapply(models, lmBF, data = data, whichRandom = whichRandom,
                    rscaleFixed = rscaleFixed, rscaleRandom = rscaleRandom,
                    iterations = iterations)
      
    }
    
    # combine all the Bayes factors into one BFBayesFactor object
    bfObj = do.call("c", bfs)
    # If we have random effects, make those the denominator
    if(length(whichRandom) > 0) bfObj = bfObj[-length(bfObj)] / bfObj[length(bfObj)]
    
    return(bfObj)
  }


