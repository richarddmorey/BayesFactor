

nWayAOV.MC<- function(y, X=NULL, struc=NULL, dataFixed=NULL, dataRandom=NULL,
           modelFormula = NULL, iterations = 10000, rscale = 1, rscaleFixed="medium", rscaleRandom=1, progress = FALSE, 
           samples = FALSE, gsamples = FALSE, gibi = NULL, logbf = FALSE)
{
  if(!is.null(X) & !is.null(struc)){
    builtDesign = FALSE
    message("Using X and struc arguments to build model.")
    if(is.null(rscale)){
      stop("rscale must be specified if using design matrix.")
    }
  }else if(!is.null(dataFixed) | !is.null(dataRandom)){
    builtDesign = TRUE
    
    rscaleFixed = rpriorValues("allNways","fixed",rscaleFixed)
    rscaleRandom = rpriorValues("allNways","random",rscaleRandom)
    
    # Convert to factors if needed.
    if(!is.null(dataFixed)){
      dataFixed = data.frame(dataFixed)
      if(any(!sapply(dataFixed,is.factor)) & !is.null(dataFixed)){
        dataFixed <- data.frame(lapply(dataFixed, as.factor))
        message("Converted columns of dataFixed to factors.")
      }
    }
    if(!is.null(dataRandom)){
      dataRandom = data.frame(dataRandom)
      if(any(!sapply(dataRandom,is.factor)) & !is.null(dataRandom)){
        dataRandom <- data.frame(lapply(dataRandom, as.factor))
        message("Converted columns of dataRandom to factors.")
      } 
    }
    ## Build design matrix of full model
    bfEnv = new.env(parent = baseenv())

    bfEnv$dataFixed = dataFixed
    bfEnv$dataRandom = dataRandom    
    
    if(is.null(modelFormula)){
      message("Using dataFixed/dataRandom to build full model.")
      modelFormula = topModelFormula(bfEnv)
    }else{
      message("Using formula and dataFixed/dataRandom to build model.")
    } 
    modInfo = buildModelInfoFormula(modelFormula, bfEnv)
    X = modInfo$X
    struc = modInfo$struc
    effType = modInfo$type

    rscale = sapply(effType, function(type){
      if(type=="random") return(rscaleRandom)
      if(type=="fixed") return(rscaleFixed)
    })
    
  }else{
    stop("Insufficient information specified to build model. Specify either X and struc, or dataFixed/dataRandom as appropriate.")
  }
  
  
  y = as.numeric(y)  
  X = as.numeric(X)
	struc = unlist(struc)
	
	X = matrix(X,nrow=length(y))
	
	
	N = as.integer(dim(X)[1])
	if(all(X[,1]==1))
	{
		P = as.integer(dim(X)[2])-1
		X = as.matrix(X[,-1],N,P)
	}else{
		P = as.integer(dim(X)[2])
	}
	
	if(sum(struc) != P)
	{
		stop(paste("Invalid struc argument. sum(struc) must be the the same as the number of parameters (excluding intercept):",sum(struc),"!=",P))
	}
	nGs = length(struc)
	
	nullLike = - ((N-1)/2)*log((N-1)*var(y))
	
	iterations = as.integer(iterations)

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
	
	C = diag(N) - matrix(1/N,N,N)
	XtCX = t(X) %*% C %*% X
	XtCy = t(X) %*% C %*% as.matrix(y,cols=1)
	ytCy = var(y)*(N-1)
	
	a = rep(0.5,nGs)
	if(length(rscale)==nGs){
		b = rscale^2/2
	}else if(length(rscale)==1){		
		b = rep(rscale^2/2,nGs)
	}else{
		stop(paste("Length of rscale vector wrong. Was",length(rscale),"and should be",nGs,"."))
	}
	
	gMap = as.integer(inverse.rle(list(values = (1:nGs)-1, lengths = struc)))
	
	returnList = .Call("RjeffSamplerNwayAov", iterations, XtCX, XtCy, ytCy, N, P, nGs, gMap, a, b,
				as.integer(iterations/100*progress), pbFun, new.env(), package="BayesFactor")

	if(inherits(pb,"txtProgressBar")) close(pb);
	
	bf = returnList[[1]] - nullLike
	if(!logbf)
	{
		bf = exp(bf)
	}
	
	if(samples & gsamples)
	{
		return(list(bf,returnList[[2]],returnList[[3]]))
	}else if(!samples & gsamples){
		return(list(bf,returnList[[3]]))
	}else if(samples & !gsamples){
		return(list(bf,returnList[[2]]))
	}else{
		return(bf)
	}
}

topModelFormula <- function(env)
{  
  dataFixed = env$dataFixed
  dataRandom = env$dataRandom
  
  if(!is.null(dataFixed)){
    fixedPart = paste(colnames(dataFixed),collapse=" * ")
  }else{
    fixedPart = NULL
  }
  if(!is.null(dataRandom)){
    randomPart = paste(colnames(dataRandom),collapse=" + ")
  }else{
    randomPart = NULL
  }
  wholeFormula = paste(c(fixedPart,randomPart),collapse=' + ')
  return(as.formula(paste("~",wholeFormula)))
}

buildModelInfoFormula <- function(formula, env)
{
  terms = attr(terms(formula),"term.labels")
  columns = unique(unlist(strsplit(terms,':',fixed=TRUE)))
  
  
  dataFixed = env$dataFixed
  dataRandom = env$dataRandom
  
  randomEffs = columns[columns %in% colnames(dataRandom)]
  fixedEffs = columns[columns %in% colnames(dataFixed)]
    
  if(any(randomEffs %in% fixedEffs)) stop("Duplicate names in dataFixed and dataRandom.")
  
  effTypes <- sapply(terms, fixedOrRandom, fixedEffs=fixedEffs, randomEffs=randomEffs) 
 
  
  randomDesigns = lapply(dataRandom,dmat1,fixed=FALSE)
  fixedDesigns = lapply(dataFixed,dmat1,fixed=TRUE)
  names(randomDesigns) = colnames(dataRandom)
  names(fixedDesigns) = colnames(dataFixed)
  allDesigns = c(fixedDesigns,randomDesigns)
  
  modelDesignList <- lapply(terms,dmatGeneric,allDesigns=allDesigns)
  struc = sapply(modelDesignList,ncol)
  X = do.call(cbind,modelDesignList)
  
  return(list(X=X,struc=struc,types=effTypes))
}

fixedOrRandom <- function(term, fixedEffs, randomEffs){
  splterm = strsplit(term,':',fixed=TRUE)[[1]]
  if(length(splterm)==1){
    if(splterm %in% randomEffs){
      return("random")
    }else if(splterm %in% fixedEffs){
      return("fixed")
    }
  }else{
    if(any(splterm %in% randomEffs)){
      return("random")
    }else if(all(splterm %in% fixedEffs)){
      return("fixed")
    }
  }
  stop(paste("Could not determine column type for ",term,"- was it included in data?"))
}

dmatGeneric <- function(term,allDesigns){
  elements = strsplit(term,':',fixed=TRUE)[[1]]
  if(length(elements)==1){
      return(allDesigns[[term]])
  }else{
      return(design.mat.intList(allDesigns[elements]))
  }
}

dmat1 <- function(data, fixed=FALSE){
  data = data.frame(x = data)
  data$x = as.factor(data$x)
  X = model.matrix(~ x - 1, data)
  X = matrix(X,nrow=length(data$x))
  if(fixed){
    X = X %*% fixedFromRandomProjection(ncol(X))
  }
  return(X)
}




nWayAOV.Gibbs <- function(y,X=NULL,struc=NULL,dataFixed=NULL, dataRandom=NULL,modelFormula = NULL,iterations=10000,rscale = 1,rscaleFixed="medium", rscaleRandom=1, progress=TRUE, unreduce=TRUE)
{
  if(!is.null(X) & !is.null(struc)){
    builtDesign = FALSE
    message("Using X and struc arguments to build model.")
    if(is.null(rscale)){
      stop("rscale must be specified if using design matrix.")
    }
  }else if(!is.null(dataFixed) | !is.null(dataRandom)){
    builtDesign = TRUE
    
    rscaleFixed = rpriorValues("allNways","fixed",rscaleFixed)
    rscaleRandom = rpriorValues("allNways","random",rscaleRandom)
    
    # Convert to factors if needed.
    if(!is.null(dataFixed)){
      dataFixed = data.frame(dataFixed)
      if(any(!sapply(dataFixed,is.factor)) & !is.null(dataFixed)){
        dataFixed <- data.frame(lapply(dataFixed, as.factor))
        message("Converted columns of dataFixed to factors.")
      }
    }
    if(!is.null(dataRandom)){
      dataRandom = data.frame(dataRandom)
      if(any(!sapply(dataRandom,is.factor)) & !is.null(dataRandom)){
        dataRandom <- data.frame(lapply(dataRandom, as.factor))
        message("Converted columns of dataRandom to factors.")
      } 
    }
    ## Build design matrix of full model
    bfEnv = new.env(parent = baseenv())
    
    bfEnv$dataFixed = dataFixed
    bfEnv$dataRandom = dataRandom    
    
    if(is.null(modelFormula)){
      message("Using dataFixed/dataRandom to build full model.")
      modelFormula = topModelFormula(bfEnv)
    }else{
      message("Using formula and dataFixed/dataRandom to build model.")
    } 
    modInfo = buildModelInfoFormula(modelFormula, bfEnv)
    X = modInfo$X
    struc = modInfo$struc
    effType = modInfo$type
    
    rscale = sapply(effType, function(type){
      if(type=="random") return(rscaleRandom)
      if(type=="fixed") return(rscaleFixed)
    })
    
  }else{
    stop("Insufficient information specified to build model. Specify either X and struc, or dataFixed/dataRandom as appropriate.")
  }
  
  
  nLevels <- apply(X,2,function(v) length(unique(v)))
  if(any(nLevels==1)){
    const = which(nLevels==1)
    if(length(const)==1){
      Z = cbind(1,X[,-const])
      #warning("Constant column in design matrix stripped.")
    }else{
      stop("Error: More than one constant column found in design matrix!")
    }
  }else{
    Z = cbind(1,X)
  }
  p = ncol(Z)-1
  N = length(y)
  if(N != nrow(Z)) stop("Error: Data and design matrix do not conform.")
  
  if (sum(struc) != p) {
    stop(paste("Invalid struc argument. sum(struc) must be the the same as the number of parameters (excluding intercept):", 
               sum(struc), "!=", p))
  }
  
  if (progress) {
    pb = txtProgressBar(min = 0, max = 100, style = 3)
  }
  else {
    pb = NULL
  }
  pbFun = function(samps){ 
    	if(progress){
    		percent = as.integer(round(samps / iterations * 100))
    	    setTxtProgressBar(pb, percent)
    	}
    }
  
  if(length(rscale)==1){
    rscale = struc*0+rscale
  }else if(length(rscale)!=length(struc)){
    stop("Error: invalid scale vector rscale.")
  }
  
  ZtZ = t(Z)%*%Z
  Zty = t(Z)%*%matrix(y,ncol=1)
  gMap = inverse.rle(list(values=1:length(struc),lengths=struc))-1


  chains = .Call("RGibbsNwayAov", as.integer(iterations), y, Z, ZtZ, Zty, as.integer(N), 
  								  as.integer(p), length(struc), as.integer(gMap), rscale,
								  as.integer(progress), pbFun, new.env(), 
								  package="BayesFactor")

  if(progress) close(pb);
	dim(chains) <- c(2 + p + length(struc), iterations)
	chains = mcmc(t(chains))  
	if(builtDesign){
    chains = makeChainNeater(chains,struc,modelFormula,dataFixed,dataRandom,unreduce=unreduce)	  
  }else{
    labels = c("mu",paste("beta",1:p,sep="_"),"sig2",paste("g",1:length(struc),sep="_"))
	  colnames(chains) = labels
  }
	return(chains)
}


makeChainNeater <- function(chains,struc,formula,dataFixed,dataRandom,unreduce){
  labels = c("mu")
  
  terms = attr(terms(formula),"term.labels")
  types = c()
  types[colnames(dataRandom)] = "random"
  if(unreduce){
    types[colnames(dataFixed)] = "random"  
  }else{
    types[colnames(dataFixed)] = "fixed"
  }
  
  if(is.null(dataFixed)){
    allData = dataRandom
  }else if(is.null(dataRandom)){
    allData = dataFixed
  }else if(!is.null(dataRandom) & !is.null(dataFixed)){
    allData = cbind(dataFixed,dataRandom)
  }else{
    stop("No data.")
  }
  
  labelList = lapply(terms, 
         function(term, types){
           terms = strsplit(term,':',fixed=TRUE)[[1]]
           my.names = design.names.intList(allData[terms], types)
           return(paste(term,"-",my.names))
          },
         types=types)
  types[colnames(dataFixed)] = "fixed"
  labels = c(labels,unlist(labelList))
  if(unreduce & any(types=="fixed")){
    betaChains = chains[,1:sum(struc) + 1]
    parMap = inverse.rle(list(values=1:length(struc),lengths=struc))
    unreducedChains = lapply(1:length(struc), 
           function(el, chains, parMap, terms, types, data){
             chains = chains[,parMap==el]
             term = strsplit(terms[el],':',fixed=TRUE)[[1]]
             S = design.projection.intList(data[term], types)
             return(chains%*%t(S))
           }, 
           chains = betaChains, parMap = parMap,
           terms = terms, types = types,
           data = allData
           )
    betaChains = do.call(cbind, unreducedChains)
    chains = cbind(chains[,1],betaChains,chains[,-(1:(sum(struc) + 1))])
  }
  
  labels = c(labels, "sig2",paste("g",1:length(struc),sep="_"))
  
  colnames(chains) = labels  
  return(chains)
}
