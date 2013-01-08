createRscales <- function(formula, data, dataTypes, rscaleFixed = NULL, rscaleRandom = NULL, rscaleCont = NULL){
  
  rscaleFixed = rpriorValues("allNways","fixed",rscaleFixed)
  rscaleRandom = rpriorValues("allNways","random",rscaleRandom)
  rscaleCont = rpriorValues("regression",,rscaleCont)
  
  types = termTypes(formula, data, dataTypes)
  nFac = sum( (types=="random") | (types=="fixed") ) 
  nCont = any(types=="continuous") * 1
  nGs = nFac + nCont
  
  rscale = 1:nGs * NA
  rscaleTypes = rscale
  
  if(nCont > 0) rscaleTypes[nGs] = "continuous"
  if(nFac > 0){
    facTypes = types[types != "continuous"]
    rscaleTypes[1:nFac] = facTypes 
  }
  
  rscale[rscaleTypes=="continuous"] = rscaleCont
  rscale[rscaleTypes=="fixed"] = rscaleFixed
  rscale[rscaleTypes=="random"] = rscaleRandom
  
  return(rscale)
}


createGMap <- function(formula, data, dataTypes){
  
  factors = fmlaFactors(formula, data)[-1]
  if(length(factors)<1) return(c())
  
  # Compute number of parameters for each specified column
  nXcols = numColsForFactor(formula, data, dataTypes)
  lvls = termLevels(formula, data, nXcols)
  types = termTypes(formula, data, dataTypes)
  
  # each random or fixed group gets a parameter, and all continuous together get 1
  nFac = sum( (types=="random") | (types=="fixed") ) 
  nCont = any(types=="continuous") * 1
  nGs = nFac + nCont
  P = sum(lvls)
  
  gGroups = inverse.rle(list(lengths=lvls,values=names(lvls)))
  gMap = 1:P * NA
  names(gMap) = gGroups
  
  gGroupsFac = lvls[types != "continuous"] * 0 + (1:nFac - 1)
  
  gMap[types[gGroups] == "continuous"] = nGs - 1
  gMap[types[gGroups] != "continuous"] = gGroupsFac[names(gMap[types[gGroups] != "continuous"])]
  
  return(gMap)
}

numColsForFactor <- function(formula, data, dataTypes){
  factors = fmlaFactors(formula, data)[-1]
  sapply(factors, function(el, data, dataTypes){
    switch(dataTypes[el],
           fixed = nlevels(data[,el]) - 1,
           random = nlevels(data[,el]),
           continuous = 1
    )
  }, data = data, dataTypes = dataTypes)
}

termLevels <- function(formula, data, nXcols){
  trms = attr(terms(formula, data = data),"term.labels")
  sapply(trms, function(term, nXcols){
    constit = strsplit(term, ":", fixed=TRUE)[[1]]
    prod(nXcols[constit])
  }, nXcols = nXcols)
}

termTypes <- function(formula, data, dataTypes){
  trms = attr(terms(formula, data = data),"term.labels")
  sapply(trms, function(term, dataTypes){
    constit = strsplit(term, ":", fixed=TRUE)[[1]]
    types = dataTypes[constit]
    if(any(types=="continuous")) return("continuous")
    if(any(types=="random")) return("random")
    return("fixed")
  }, dataTypes = dataTypes)
}


fullDesignMatrix <- function(fmla, data, dataTypes){
  trms <- attr(terms(fmla, data = data), "term.labels")
  
  Xs = lapply(trms,function(trm, data, dataTypes){
    oneDesignMatrix(trm, data = data, dataTypes = dataTypes)      
  }, data = data, dataTypes = dataTypes)    
  
  do.call("cbind",Xs)
}

oneDesignMatrix <- function(trm, data, dataTypes)
{
  effects <- unlist(strsplit(trm, ":", fixed = TRUE))
  #check to ensure all terms are in data
  checkEffects(effects, data, dataTypes)
  
  if(length(effects) == 1){
    effect = paste("~",effects,"-1")
    X = model.matrix(formula(effect),data = data)
    if(dataTypes[effects] == "fixed"){
      X = X %*% fixedFromRandomProjection(ncol(X))
      colnames(X) = paste(effects,"_redu_",1:ncol(X),sep="")
    }
    return(X)
  }else{
    Xs = lapply(effects, function(trm, data, dataTypes){
      oneDesignMatrix(trm, data = data, dataTypes = dataTypes)
    }, data = data, dataTypes = dataTypes)
    X = Reduce(rowMultiply, x = Xs)
    return(X)
  }
}

design.names.intList <- function(effects, data, dataTypes){
  type = dataTypes[ effects[1] ]
  firstCol = data[ ,effects[1] ]
  nLevs = nlevels( firstCol )
  if(length(effects)==1){
    if(type=="random") return(levels(firstCol))
    if(type=="fixed") return(0:(nLevs-2))
    if(type=="continuous") return(effects)
  }else{
    if(type=="random") 
      return(rowPaste(levels(firstCol), design.names.intList(effects[-1], data, dataTypes) ))
    if(type=="fixed") 
      return(rowPaste(0:(nLevs-2), design.names.intList(effects[-1], data, dataTypes) ))
    if(type=="continuous") 
      return(rowPaste(0:(nLevs-2), design.names.intList(effects[-1], data, dataTypes) ))
  }    
}

design.projection.intList <- function(effects, data, dataTypes){
  type = dataTypes[ effects[1] ]
  firstCol = data[ ,effects[1] ]
  nLevs = nlevels( firstCol )
  if(length(effects)==1){
    if(type=="random" | type=="continuous") return(diag(nLevs))
    if(type=="fixed") return(fixedFromRandomProjection(nLevs))
  }else{
    if(type=="random") 
      return(kronecker(diag(nLevs), design.projection.intList(effects[-1],data, dataTypes) ))
    if(type=="fixed") 
      return(kronecker(fixedFromRandomProjection(nLevs), design.projection.intList(effects[-1], data, dataTypes) ))
    if(type=="continuous") 
      return( design.projection.intList(effects[-1], data, dataTypes) )
  }    
}

rowPaste = function(v1,v2)
{
  as.vector(t(outer(v1,v2,paste,sep=".&.")))
}

rowMultiply = function(x,y)
{
  if(nrow(x) != nrow(y)) stop("Unequal row numbers in row.multiply:", nrow(x),", ",nrow(y))
  K = sapply(1:nrow(x), function(n, x, y){
    kronecker(x[n,], y[n,])
  }, x = x, y = y )
  # add names
  K <- t(matrix(K, ncol = nrow(x)))
  colnames(K) = as.vector(t(
    outer(colnames(x), colnames(y), function(x,y){
      paste(x, y,sep=".&.")
    })))
  return(K)
}

# Create projection matrix
fixedFromRandomProjection <- function(nlevRandom){
  centering=diag(nlevRandom)-(1/nlevRandom)
  S=(eigen(centering)$vectors)[,1:(nlevRandom-1)]
  return(matrix(S,nrow=nlevRandom))
}


nWayFormula <- function(formula, data, dataTypes, rscaleFixed=NULL, rscaleRandom=NULL, rscaleCont=NULL, gibbs=FALSE, unreduce=TRUE, ...){
  
  checkFormula(formula, data, analysis = "lm")
  y = data[,deparse(formula[[2]])]
  
  X = fullDesignMatrix(formula, data, dataTypes)
  
  rscale = createRscales(formula, data, dataTypes, rscaleFixed, rscaleRandom, rscaleCont)
  gMap = createGMap(formula, data, dataTypes)
  
  retVal = nWayAOV(y, X, gMap = gMap, rscale = rscale, gibbs = gibbs, ...)
  if(gibbs){
    retVal <- mcmc(makeChainNeater(retVal, colnames(X), formula, data, dataTypes, gMap, unreduce))  
  }
  return(retVal)
}  


makeLabelList <- function(formula, data, dataTypes, unreduce){
  
  terms = attr(terms(formula, data = data), "term.labels")
  
  if(unreduce) 
    dataTypes[dataTypes == "fixed"] = "random"
  
  labelList = lapply(terms, 
                     function(term, data, dataTypes){
                       effects = strsplit(term,":",fixed=TRUE)[[1]]
                       my.names = design.names.intList(effects, data, dataTypes)
                       return(paste(term,"-",my.names,sep=""))
                     },
                     data = data, dataTypes=dataTypes)
  
  # join them all together in one cector
  unlist(labelList)
}

unreduceChainPart = function(term, chains, data, dataTypes, gMap){
  effects = strsplit(term,":", fixed = TRUE)[[1]]
  chains = chains[, names(gMap)==term, drop = FALSE ]
  if(any(dataTypes[effects]=="fixed")){
    S = design.projection.intList(effects, data, dataTypes)
    return(chains%*%t(S))
  }else{
    return(chains)
  }
}

ureduceChains = function(chains, formula, data, dataTypes, gMap){
  
  terms = attr(terms(formula, data = data), "term.labels")
  
  unreducedChains = lapply(terms, unreduceChainPart, chains=chains, data = data, dataTypes = dataTypes, gMap = gMap)  
  do.call(cbind, unreducedChains)
}

makeChainNeater <- function(chains, Xnames, formula, data, dataTypes, gMap, unreduce){
  P = length(gMap)
  nGs = max(gMap) + 1
  factors = fmlaFactors(formula, data)[-1]
  dataTypes = dataTypes[ names(dataTypes) %in% factors ]
  
  if(!unreduce | !any(dataTypes == "fixed")) {
    labels = c("mu", Xnames, "sig2", paste("g",1:nGs,sep="_"))
    colnames(chains) = labels 
    return(chains)
  }
  
  labels = c("mu")  
  betaChains = chains[,1:P + 1, drop = FALSE]
  types = termTypes(formula, data, dataTypes)
  
  # Make column names 
  parLabels = makeLabelList(formula, data, dataTypes, unreduce)
  labels = c(labels, parLabels)
  
  betaChains = ureduceChains(betaChains, formula, data, dataTypes, gMap)
  
  newChains = cbind(chains[,1],betaChains,chains[,-(1:(P + 1))])
  
  labels = c(labels, "sig2",paste("g",1:nGs,sep="_"))
  colnames(newChains) = labels 
  return(newChains)
}
