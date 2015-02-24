singleGBayesFactor <- function(y,X,rscale,gMap){
  if(ncol(X)==1){
    dat = data.frame(y=y,x=as.factor(X[,1])) 
    freqs = table(dat$x)
    t = t.test(y~x,data=dat, var.eq=TRUE)$statistic
    bf = ttest.tstat(t=t, n1=freqs[1], n2=freqs[2],rscale=rscale*sqrt(2))
    return(bf)
  }else{
    
    incCont = 0
    priorX = matrix(1,0,0)
    N = length(y)
    Cny = matrix(y - mean(y), N)
    CnX = t(t(X) - colMeans(X))
    CnytCnX = t(Cny)%*%CnX
    sumSq = var(y) * (N-1) 
    gMapCounts = table(gMap)
    
    f1 = Vectorize(
      function(g,...){
        exp(Qg(log(g), ..., limit=FALSE) - log(g) - const)
      },"g")
    
    integral = try({
      op = optim(0, Qg, control=list(fnscale=-1),gr=dQg, method="BFGS",
                 sumSq=sumSq,Cny=Cny,CnX=CnX,CnytCnX=CnytCnX, rscale=rscale, gMap=gMap, gMapCounts=gMapCounts,priorX=priorX,incCont=incCont)
      const = op$value - op$par
      integrate(f1,0,Inf,sumSq=sumSq,Cny=Cny,CnX=CnX,CnytCnX=CnytCnX,rscale=rscale,gMap=gMap,gMapCounts=gMapCounts,const=const,priorX=priorX,incCont=incCont)
    })
    if(inherits(integral,"try-error")){
      return(list(bf = NA, properror = NA))
    }
    lbf = log(integral$value)
    prop.error = exp(log(integral$abs.error) - lbf)
    return(list(bf = lbf + const, properror = prop.error))
  }
}


doNwaySampling<-function(method, y, X, rscale, iterations, gMap, incCont, progress, callback = function(...) as.integer(0))
{
  goodSamples = NULL
  simpSamples = NULL
  impSamples = NULL
  apx = NULL
  testNsamples = options()$BFpretestIterations
  
  testCallback = function(...) callback(0)
  
  if(ncol(X)==1) method="simple"
  
  if(method=="auto"){
    simpSamples = try(jzs_sampler(testNsamples, y, X, rscale, gMap, incCont, NA, NA, FALSE, testCallback, 1, 0))
    simpleErr = propErrorEst(simpSamples)
    logAbsSimpErr = logMeanExpLogs(simpSamples) + log(simpleErr) 
     
    
    apx = suppressWarnings(try(gaussianApproxAOV(y,X,rscale,gMap,incCont)))
    if(inherits(apx,"try-error")){
      method="simple"
    }else{
      impSamples = try(jzs_sampler(testNsamples, y, X, rscale, gMap, incCont, apx$mu, apx$sig, FALSE, testCallback, 1, 1))
      if(inherits(impSamples, "try-error")){
        method="simple"
      }else{
        impErr = propErrorEst(impSamples)
        logAbsImpErr = logMeanExpLogs(impSamples) + log(impErr)
        if(is.na(impErr)){
          method="simple"
        }else if(is.na(simpleErr)){
          method="importance"
        }else{
          method = ifelse(impErr>simpleErr,"simple","importance")
        }
      }  
    }
  }
  
  if(method=="importance"){

    if(is.null(apx) | inherits(apx,"try-error"))  
      apx = try(gaussianApproxAOV(y,X,rscale,gMap,incCont))
    if(inherits(apx, "try-error")){
      method="simple"   
    }else{
      goodSamples= try(jzs_sampler(iterations, y, X, rscale, gMap, incCont, apx$mu, apx$sig, progress, callback, 1, 1))
      if(inherits(goodSamples,"try-error")){
        method="simple"
        goodSamples = NULL
      }
    }
  }  
  if(method=="simple" | is.null(goodSamples)){
    method = "simple"
    goodSamples = jzs_sampler(iterations, y, X, rscale, gMap, incCont, NA, NA, progress, callback, 1, 0)
    
  }
  if(is.null(goodSamples)){  
    warning("Unknown sampling method requested (or sampling failed) for nWayAOV")
    return(c(bf=NA,properror=NA))
  }
  
  if( any(is.na(goodSamples)) ) warning("Some NAs were removed from sampling results: ",sum(is.na(goodSamples))," in total.")
  bfSamp = goodSamples[!is.na(goodSamples)]
  n2 = length(bfSamp)
    
  bf = logMeanExpLogs(bfSamp)
  
  # estimate error
  properror = propErrorEst(bfSamp)
  
  return(list(bf = bf, properror=properror, N = n2, method = method, sampled = TRUE, code = randomString(1)))
}

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
  
  sparse = any(dataTypes=="random")
  
  Xs = lapply(trms,function(trm, data, dataTypes){
    oneDesignMatrix(trm, data = data, dataTypes = dataTypes, sparse = sparse)      
  }, data = data, dataTypes = dataTypes)    
  
  do.call("cBind" ,Xs)
}

oneDesignMatrix <- function(trm, data, dataTypes, sparse = FALSE)
{
  effects <- unlist(strsplit(trm, ":", fixed = TRUE))
  #check to ensure all terms are in data
  checkEffects(effects, data, dataTypes)
    
  if(length(effects) == 1){
    effect = paste("~",effects,"-1")
    X = model.Matrix(formula(effect),data = data, sparse = sparse)
    if(dataTypes[effects] == "fixed"){
      X = X %*% fixedFromRandomProjection(ncol(X), sparse = sparse)
      colnames(X) = paste(effects,"_redu_",1:ncol(X),sep="")
    }
    return(X)
  }else{
    Xs = lapply(effects, function(trm, data, dataTypes, sparse){
      oneDesignMatrix(trm, data = data, dataTypes = dataTypes, sparse = sparse)
    }, data = data, dataTypes = dataTypes, sparse = sparse)
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
      return( design.names.intList(effects[-1], data, dataTypes) )
      #return(rowPaste(0:(nLevs-2), design.names.intList(effects[-1], data, dataTypes) ))
  }    
}

design.projection.intList <- function(effects, data, dataTypes){
  type = dataTypes[ effects[1] ]
  firstCol = data[ ,effects[1] ]
  nLevs = nlevels( firstCol )
  if(length(effects)==1){
    if(type=="random") return(bdiag(diag(nLevs)))
    if(type=="fixed") return(fixedFromRandomProjection(nLevs))
    if(type=="continuous") return(Matrix(1,1,1))
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

rowMultiply = function(x, y)
{
  sparse = is(x, "sparseMatrix") | is(y, "sparseMatrix")
  if(nrow(x) != nrow(y)) stop("Unequal row numbers in row.multiply:", nrow(x),", ",nrow(y))
  K = sapply(1:nrow(x), function(n, x, y){
    kronecker(x[n,], y[n,])
  }, x = x, y = y )
  # add names
  K <- t(Matrix(as.vector(K), ncol = nrow(x), sparse = sparse))
  colnames(K) = as.vector(t(
    outer(colnames(x), colnames(y), function(x,y){
      paste(x, y,sep=".&.")
    })))
  return(K)
}

# Create projection matrix
fixedFromRandomProjection <- function(nlevRandom, sparse = FALSE){
  centering=diag(nlevRandom)-(1/nlevRandom)
  S=as.vector((eigen(centering)$vectors)[,1:(nlevRandom-1)])
  return(Matrix(S,nrow=nlevRandom, sparse = sparse))
}

centerContinuousColumns <- function(data){
  mycols = lapply(data,function(colmn){
    if(is.factor(colmn)){
      return(colmn)
    }else{
      return(colmn - mean(colmn))
    }
  })
  return(data.frame(mycols))
}

nWayFormula <- function(formula, data, dataTypes, rscaleFixed=NULL, rscaleRandom=NULL, rscaleCont=NULL, posterior=FALSE, columnFilter = NULL, unreduce=TRUE, ...){
  
  checkFormula(formula, data, analysis = "lm")
  
  
  y = data[,stringFromFormula(formula[[2]])]
  data <- centerContinuousColumns(data)
  X = fullDesignMatrix(formula, data, dataTypes)
  
  # To be removed when sparse matrix support is complete
  X = as.matrix(X)
  
  rscale = createRscales(formula, data, dataTypes, rscaleFixed, rscaleRandom, rscaleCont)
  gMap = createGMap(formula, data, dataTypes)
  
  if(any(dataTypes=="continuous")){
    continuous = termTypes(formula, data, dataTypes)=="continuous"
    continuous = continuous[names(gMap)]
  }else{
    continuous = FALSE
  }
  
  ## Determine which columns we will ignore
  if(is.null(columnFilter)){
    ignoreCols = NULL
  }else{
    ignoreCols = filterVectorLogical(columnFilter, names(gMap))
  }
  if(all(ignoreCols) & !is.null(ignoreCols)) stop("Filtering out all chain columns of interest is not allowed.")
  
  retVal = nWayAOV(y, X, gMap = gMap, rscale = rscale, posterior = posterior, continuous = continuous, ignoreCols=ignoreCols,...)
  if(posterior){
    retVal <- mcmc(makeChainNeater(retVal, colnames(X), formula, data, dataTypes, gMap, unreduce, continuous, columnFilter))  
  }
  return(retVal)
}  


makeLabelList <- function(formula, data, dataTypes, unreduce, columnFilter){
  
  terms = attr(terms(formula, data = data), "term.labels")
  if(!is.null(columnFilter))
    terms = terms[!filterVectorLogical(columnFilter,terms)]
  
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

unreduceChainPart = function(term, chains, data, dataTypes, gMap, ignoreCols){
  effects = strsplit(term,":", fixed = TRUE)[[1]]
  myCols = names(gMap)==term
  if(ignoreCols[myCols][1]) return(NULL)
  
  # Figure out which columns we need to look at, given that some are missing
  cumulativeIgnored = sum(ignoreCols[1:which(myCols)[1]]) # How many are ignored up to the one of interest?
  remappedCols = which(myCols) - cumulativeIgnored
  chains = chains[, remappedCols, drop = FALSE ]
  if(any(dataTypes[effects]=="fixed")){
    S = design.projection.intList(effects, data, dataTypes)
    return(chains%*%as.matrix(t(S)))
  }else{
    return(chains)
  }
}

ureduceChains = function(chains, formula, data, dataTypes, gMap, ignoreCols){
  
  terms = attr(terms(formula, data = data), "term.labels")
  
  unreducedChains = lapply(terms, unreduceChainPart, chains=chains, data = data, dataTypes = dataTypes, gMap = gMap, ignoreCols=ignoreCols)  
  do.call(cbind, unreducedChains)
}

makeChainNeater <- function(chains, Xnames, formula, data, dataTypes, gMap, unreduce, continuous, columnFilter){
  P = length(gMap)
  nGs = max(gMap) + 1
  factors = fmlaFactors(formula, data)[-1]
  dataTypes = dataTypes[ names(dataTypes) %in% factors ]
  types = termTypes(formula, data, dataTypes)
   
  lastPars = ncol(chains) + (-nGs):0
  
  if(any(continuous)){
    gNames = paste("g",c(names(types[types!="continuous"]),"continuous"),sep="_")
  }else{
    gNames = paste("g",names(types), sep="_")  
  }
    
  if(is.null(columnFilter)){
    ignoreCols = ignoreCols = rep(0,P)
  }else{
    ignoreCols=filterVectorLogical(columnFilter, names(gMap))
  }
  
  if(!unreduce | !any(dataTypes == "fixed")) {
    labels = c("mu", Xnames[!ignoreCols], "sig2", gNames)
    colnames(chains) = labels 
    return(chains)
  }
  
  # Make column names 
  parLabels = makeLabelList(formula, data, dataTypes, unreduce, columnFilter)
    
  labels = c("mu", parLabels)

  betaChains = chains[,1:(ncol(chains)-2-nGs) + 1, drop = FALSE]  
  betaChains = ureduceChains(betaChains, formula, data, dataTypes, gMap, ignoreCols)
  
  newChains = cbind(chains[,1],betaChains,chains[,lastPars])
  
  labels = c(labels, "sig2", gNames)
  colnames(newChains) = labels 
  return(newChains)
}

