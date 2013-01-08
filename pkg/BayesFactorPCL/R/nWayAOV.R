

createRscales <- function(formula, dataTypes, rscaleFixed = NULL, rscaleRandom = NULL, rscaleCont = NULL){
  
  rscaleFixed = rpriorValues("allNways","fixed",rscaleFixed)
  rscaleRandom = rpriorValues("allNways","random",rscaleRandom)
  rscaleCont = rpriorValues("regression",,rscaleCont)
  
  types = termTypes(formula, dataTypes)
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
  
  factors = fmlaFactors(formula)[-1]
  if(length(factors)<1) return(c())
  
  # Compute number of parameters for each specified column
  nXcols = numColsForFactor(formula, data, dataTypes)
  lvls = termLevels(formula, nXcols)
  types = termTypes(formula, dataTypes)
  
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
  factors = fmlaFactors(formula)[-1]
  sapply(factors, function(el, data, dataTypes){
    switch(dataTypes[el],
           fixed = nlevels(data[,el]) - 1,
           random = nlevels(data[,el]),
           continuous = 1
    )
  }, data = data, dataTypes = dataTypes)
}

termLevels <- function(formula, nXcols){
  trms = attr(terms(formula),"term.labels")
  sapply(trms, function(term, nXcols){
    constit = strsplit(term, ":", fixed=TRUE)[[1]]
    prod(nXcols[constit])
  }, nXcols = nXcols)
}

termTypes <- function(formula, dataTypes){
  trms = attr(terms(formula),"term.labels")
  sapply(trms, function(term, dataTypes){
    constit = strsplit(term, ":", fixed=TRUE)[[1]]
    types = dataTypes[constit]
    if(any(types=="continuous")) return("continuous")
    if(any(types=="random")) return("random")
    return("fixed")
  }, dataTypes = dataTypes)
}


fullDesignMatrix <- function(fmla, data, dataTypes){
  trms <- attr(terms(fmla), "term.labels")
  
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
  
  rscale = createRscales(formula, dataTypes, rscaleFixed, rscaleRandom, rscaleCont)
  gMap = createGMap(formula, data, dataTypes)
  
  retVal = nWayAOV(y, X, gMap = gMap, rscale = rscale, gibbs = gibbs, ...)
  if(gibbs){
    retVal <- mcmc(makeChainNeater(retVal, colnames(X), formula, data, dataTypes, gMap, unreduce))  
  }
  return(retVal)
}  

nWayAOV<- function(y, X, struc = NULL, gMap = NULL, rscale, iterations = 10000, progress = FALSE, gibi = NULL, gibbs = FALSE)
{
  if(!is.numeric(y)) stop("y must be numeric.")  
  if(!is.numeric(X)) stop("X must be numeric.")  
  
  N = length(y)
	X = matrix( X, nrow=N )
  
  constantCols = apply(X,2,function(v) length(unique(v)))==1
  if( sum(constantCols) > 0 )
	{
    X = X[,-which(constantCols)]
	  warning( sum(constantCols)," constant columns removed from X." )
  }
	P = ncol(X)
  
  if(!is.null(gMap)){
    if(length(gMap) != P)
      stop("Invalid gMap argument. length(gMap) must be the the same as the number of parameters (excluding intercept): ",sum(gMap)," != ",P)
    if( !all(0:max(gMap) %in% unique(gMap)) )
      stop("Invalid gMap argument: no index can be skipped.")
    
    nGs = as.integer(max(gMap) + 1)
  }else if(!is.null(struc)){
    struc = unlist(struc)
    if(sum(struc) != P)
      stop("Invalid struc argument. sum(struc) must be the the same as the number of parameters (excluding intercept): ",sum(struc)," != ",P)
    
    nGs = length(struc)
    gMap = as.integer(inverse.rle(list(values = (1:nGs)-1, lengths = struc)))
  }else{
    stop("One of gMap or struc must be defined.")
  }
  
  C = diag(N) - matrix(1/N,N,N)
  XtCX = t(X) %*% C %*% X
  XtCy = t(X) %*% C %*% as.matrix(y,cols=1)
  ytCy = var(y)*(N-1)
  
  a = rep(0.5,nGs)
  if(length(rscale)==nGs){
    b = rscale^2/2
  }else{
    stop("Length of rscale vector wrong. Was ", length(rscale), " and should be ", nGs,".")
  }
  
	nullLike = - ((N-1)/2)*log((N-1)*var(y))
	
  
  ####### Progress bar stuff
	if(!is.null(gibi)) {
		progress=TRUE;
		if(!is.function(gibi))
			stop("Malformed GIBI argument (not a function). You should not set this argument if running oneWayAOV.Gibbs from the console.")
	}
	if(progress & is.null(gibi)){
		pb = txtProgressBar(min = 0, max = 100, style = 3) 
	}else{ 
		pb=NULL 
	}
	
  pbFun = function(samps){ 
    if(progress){
    	percent = as.integer(round(samps / iterations * 100))
    	if(is.null(gibi)){
    		setTxtProgressBar(pb, percent)
    	}else{
    		gibi(percent)
    	}
    }
  }
  ############ End progress bar stuff
  
  
  if(gibbs){ # Create chains
    Z = cbind(1,X)
    ZtZ = t(Z)%*%Z
    Zty = t(Z)%*%matrix(y,ncol=1)
    
    chains = .Call("RGibbsNwayAov", 
                   as.integer(iterations), y, Z, ZtZ, Zty, as.integer(N), 
                   as.integer(P), as.integer(nGs), as.integer(gMap), rscale,
                   as.integer(iterations/100*progress), pbFun, new.env(), 
                   package="BayesFactor")
    
    dim(chains) <- c(2 + P + nGs, iterations)
    chains = mcmc(t(chains))  
    labels = c("mu",paste("beta",1:P,sep="_"),"sig2",paste("g",1:nGs,sep="_"))
    colnames(chains) = labels
    retVal = chains
 
  }else{ # Compute Bayes factor
    returnList = .Call("RjeffSamplerNwayAov", 
                     as.integer(iterations), XtCX, XtCy, ytCy, 
                     as.integer(N), 
                     as.integer(P), 
                     as.integer(nGs), 
                     as.integer(gMap), a, b,
                     as.integer(iterations/100*progress), 
                     pbFun, new.env(), 
                     package="BayesFactor")

    bf = returnList[[1]] - nullLike
    
    # estimate error
    bfSamp = returnList[[2]] - nullLike
    sumXsq = exp(logMeanExpLogs(2 * bfSamp)+log(iterations))
    SSq = sumXsq - exp(2*(returnList[[1]] - nullLike) + log(iterations))
    err = .5 * log(SSq) - log(iterations)
    properror = exp(err - bf)
    
    retVal = c(bf = bf, properror=properror)
  }
  
  if(inherits(pb,"txtProgressBar")) close(pb);
  return(retVal)
}

makeLabelList <- function(formula, data, dataTypes, unreduce){
  
  terms = attr(terms(formula), "term.labels")
  
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
  
  terms = attr(terms(formula), "term.labels")
  
  unreducedChains = lapply(terms, unreduceChainPart, chains=chains, data = data, dataTypes = dataTypes, gMap = gMap)  
  do.call(cbind, unreducedChains)
}

makeChainNeater <- function(chains, Xnames, formula, data, dataTypes, gMap, unreduce){
  P = length(gMap)
  nGs = max(gMap) + 1
  factors = fmlaFactors(formula)[-1]
  dataTypes = dataTypes[ names(dataTypes) %in% factors ]
    
  if(!unreduce | !any(dataTypes == "fixed")) {
    labels = c("mu", Xnames, "sig2", paste("g",1:nGs,sep="_"))
    colnames(chains) = labels 
    return(chains)
  }

  labels = c("mu")  
  betaChains = chains[,1:P + 1, drop = FALSE]
  types = termTypes(formula, dataTypes)
  
  # Make column names 
  parLabels = makeLabelList(formula, data, dataTypes, unreduce)
  labels = c(labels, parLabels)
  
  betaChains = ureduceChains(betaChains, formula, data, dataTypes, gMap)
  
  newChains = cbind(chains[,1],betaChains,chains[,-(1:(P + 1))])
  
  labels = c(labels, "sig2",paste("g",1:nGs,sep="_"))
  colnames(newChains) = labels 
  return(newChains)
}
