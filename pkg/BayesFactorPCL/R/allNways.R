allNways = function(y, dataFixed=NULL, dataRandom=NULL, iterations = 10000, 
                    whichModels="withmain", progress=TRUE, rscaleFixed="medium", 
                    rscaleRandom=1, logbf=FALSE, multicore=FALSE, extraInfo=FALSE, 
                    only.top=NULL, ... )
{
  if(!is.null(only.top)){
    warning("only.top has been deprecated; use whichModels instead. See the ?allNways for details.")
  }
  if( !(whichModels %in% c("all","top","withmain"))){
    stop("Invalid value for whichModels: must be 'all', 'top', or 'withmain'")
  }
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
  
  rscaleFixed = rpriorValues("allNways","fixed",rscaleFixed)
  rscaleRandom = rpriorValues("allNways","random",rscaleRandom)
  
  nFac = dim(dataFixed)[2]
  if(nFac==1) whichModels='all'
  bfEnv = new.env(parent = baseenv())
  designs = list()
  designs[[2^nFac]] = matrix(nrow=0,ncol=0)

  bfEnv$designMatrices = designs
  bfEnv$dataFixed = dataFixed
  bfEnv$y = y
  bfEnv$totalN = length(as.vector(y))
  bfEnv$dataRandom = dataRandom
  bfEnv$nFac = nFac

  
  if(multicore){
    allResults <- all.Nways.env.mc(env=bfEnv,iterations=iterations, whichModels, progress=FALSE, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom, ...)
  }else{
    allResults <- all.Nways.env(env=bfEnv,iterations=iterations, whichModels, progress=progress, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom, ...)
  }

  bfs = data.frame(t(allResults), stringsAsFactors=FALSE)
  

  
  topModel = ((2^(2^nFac-1))-1)
  topModelName= paste(joined.design(topModel,env=bfEnv,other="n"),collapse=" + ")
  topModelIndex = which(bfs$model==topModelName)
  
  bfs = rbind(bfs,c(bf=0, model="null", number=0,
              nParFixed=0,nParRandom=bfs$nParRandom[1],
              rscaleFixed=NA,
              rscaleRandom=bfs$rscaleRandom[1],
              omitted=topModelName))
  rownames(bfs) = bfs$model
  
  bfs$bf = as.numeric(bfs$bf)
  bfs$number = as.integer(bfs$number)
  bfs$nParFixed = as.integer(bfs$nParFixed)
  bfs$nParRandom = as.integer(bfs$nParRandom)  
  bfs$rscaleFixed = as.numeric(bfs$rscaleFixed)
  bfs$rscaleRandom = as.numeric(bfs$rscaleRandom)  
  	
  if(!is.null(dataRandom))
  {
  	nullMod = as.numeric(
                nWayAOV2(0,bfEnv,iterations=iterations, whichModels, 
                      rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom, ...)[1]
                )
  	bfs$bf = bfs$bf - nullMod
	bfs$bf[bfs$model=="null"] = 0
  }
  
  if(whichModels=="top"){
  	bfs$bf = bfs$bf - bfs$bf[topModelIndex]
  }
  
  # Prepare for return
  # sort
  bfs = bfs[order(bfs$bf, decreasing=TRUE),]
  
  # exponentiate
  if(!logbf){
  	bfs$bf = exp(bfs$bf)
  }
  
  if(extraInfo){
    return(bfs)
  }else{
    retVec = bfs$bf
    names(retVec) = bfs$model
    return(retVec)
  }
}

makeModelsVector <- function(whichModels, nFac){
  topMod = ((2^(2^nFac-1))-1)
  if(whichModels=="all"){
    modNums = 1:topMod
  }else if(whichModels=="top"){
    nDig = 2^nFac-1
    mods <- c(colSums((1-diag(nDig))*2^(0:(nDig-1))),topMod)
    modNums <- as.list(mods)
  }else if( whichModels=="withmain"){
    modNums = makeModelsAllLevels(nFac)
  }
  return(modNums)
}

all.Nways.env = function(env, whichModels, progress=progress, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom,...){
  data = env$dataFixed
	nFac = env$nFac
  modNums <- makeModelsVector(whichModels,nFac)
  
  if(progress){
    results <- pbsapply(modNums,nWayAOV2,env=env, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom,...)     
  }else{
    results <- sapply(modNums,nWayAOV2,env=env, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom,...)
  }
	return(results)
}


#### Multi core version
all.Nways.env.mc = function(env, whichModels,progress=progress, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom,...){
  if(!require(doMC)){
    warning("Required package (doMC) missing for multicore functionality. Falling back to single core functionality.")
    allResults <- all.Nways.env(env=env, whichModels=whichModels, progress=progress, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom, ...)
    return(allResults)
  }  
  
  registerDoMC()
  if(getDoParWorkers()==1){
    warning("Multicore specified, but only using 1 core. Set options(cores) to something >1.")
  }
  
  data = env$dataFixed
  nFac = env$nFac
  
  modNums <- makeModelsVector(whichModels,nFac)
  
  # No progress bar
  results <- foreach(gIndex=modNums,.combine='cbind', .options.multicore=mcoptions) %dopar% nWayAOV2(gIndex,env=env, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom,...) 
  
  return(results)
}

buildModelInfo <- function(modNum, env) {
  nFac = env$nFac
  nModels = 2^(2^nFac - 1)
  X = joined.design(modNum,env=env)
  my.name = paste(joined.design(modNum, env=env, other="n"), collapse=" + ")
  my.omitted = paste(joined.design(nModels - modNum - 1, env=env, other="n"), collapse=" + ")
  g.groups = unlist(joined.design(modNum,env=env,other="g"))
  dataRandom = env$dataRandom
  if(!is.null(dataRandom)){
    gr.groups = unlist(lapply(data.frame(dataRandom),nlevels))
  }else{
    gr.groups = NULL
  } 
  return(list(
    X = X,
    struc = c(g.groups,gr.groups),
    g.groups = g.groups,
    gr.groups = gr.groups,
    names=my.name,
    omitted=my.omitted
    ))
}


nWayAOV2 = function(modNum,env, rscaleFixed, rscaleRandom, progressCallback=NULL, ...)
{
  modInfo = buildModelInfo(modNum, env)
  y = env$y
  
  aovResults = suppressMessages(nWayAOV.MC(y, X=modInfo$X, struc=modInfo$struc, samples=FALSE,logbf=TRUE, progress=FALSE, 
                              rscale = c(rscaleFixed + modInfo$g.groups*0, rscaleRandom=rscaleRandom + modInfo$gr.groups*0),...))
  
  bfs = c(aovResults, 		
          modInfo$names,
          modNum,
          sum(modInfo$g.groups),
          sum(modInfo$gr.groups),
          ifelse(length(rscaleFixed)==1,rscaleFixed,NA),
          ifelse(length(rscaleRandom)==1,rscaleRandom,NA),
          modInfo$omitted
          )
  names(bfs) = c("bf","model","number","nParFixed","nParRandom","rscaleFixed","rscaleRandom","omitted")
  
  if(is.function(progressCallback)){
    progressCallback()
  }
  
  return(bfs)
}

unreduceChains = function(g.groups, env, chains){
  iterations <- dim(chains)[1]
  fixed.chains = chains[,2:(1+sum(g.groups))]
  fixed.chains = matrix(fixed.chains,nrow=iterations)

  stdChains = lapply(1:length(g.groups), 
    function(el, g.groups, env, chains){
      sums = c(0,cumsum(g.groups))
      C = chains[ ,sums[el] + 1:(g.groups[el]) ]
      S = other.design(el,fixed=TRUE,env,type="c")
      C%*%t(S)
    }, g.groups=g.groups, env=env, chains=fixed.chains)
  unreduced = data.frame(stdChains)
  mcmc(data.frame(chains[,1],unreduced,chains[,-(1:(1+sum(g.groups)))]))
}

design.mat.single=function(v,reduce=TRUE)
{
  nlev = nlevels(as.factor(v))
  v = as.integer(as.factor(v))
  X = matrix(0,nrow=length(v),ncol=nlev)
  X[1:length(v) + length(v)*(v-1)] = 1
  if(reduce){
    return(X %*% fixedFromRandomProjection(nlev) )
  }else{
    return(X)
  }
}

# Create the design matrix for an interaction, given the design matrices for the two factors
design.mat.int = function(X1,X2)
{
  p1 = dim(X1)[2]
  p2 = dim(X2)[2]
  X3 = apply(cbind(X1,X2),1,do.row.mult,p1=p1,p2=p2)  
  return(t(matrix(X3,ncol=dim(X1)[1])))
}

design.mat.intList <- function(Xlist){
  if(length(Xlist)==1){
    return(Xlist[[1]])
  }else{
    return(design.mat.int(Xlist[[1]], design.mat.intList(Xlist[-1]) ))
  }    
}

design.names.intList <- function(Dlist,types){
  type = types[colnames(Dlist)[1]]
  nLevs = nlevels(Dlist[[1]])
  if(length(Dlist)==1){
    if(type=="random") return(levels(Dlist[[1]]))
    if(type=="fixed") return(0:(nLevs-2))
  }else{
    if(type=="random") 
      return(do.row.paste(levels(Dlist[[1]]), design.names.intList(Dlist[-1], types) ))
    if(type=="fixed") 
      return(do.row.paste(0:(nLevs-2), design.names.intList(Dlist[-1], types) ))
  }    
}

design.projection.intList <- function(Dlist,types){
  type = types[colnames(Dlist)[1]]
  nLevs = nlevels(Dlist[[1]])
  if(length(Dlist)==1){
    if(type=="random") return(diag(nLevs))
    if(type=="fixed") return(fixedFromRandomProjection(nLevs))
  }else{
    if(type=="random") 
      return(kronecker(diag(nLevs), design.projection.intList(Dlist[-1], types) ))
    if(type=="fixed") 
      return(kronecker(fixedFromRandomProjection(nLevs), design.projection.intList(Dlist[-1], types) ))
  }    
}

do.row.paste = function(v1,v2)
{
  as.vector(t(outer(v1,v2,paste,sep=".&.")))
}

do.row.mult = function(v,p1,p2)
{
  v1 = v[1:p1]
  v2 = v[p1 + 1:p2]
  as.vector(t(outer(v1,v2,'*')))
}

joined.design = function(modelNum, env, other=NULL)
{
  N = get("totalN",env)
  dataRandom = get("dataRandom",env)
  if(modelNum==0 & is.null(dataRandom)){
      return(matrix(1,nrow=N))
  }else{
    effects = binary(modelNum)$dicotomy
    effNums = which(effects)
    if(!is.null(other)){
        gn = sapply(effNums,other.design,fixed = TRUE,env = env,type=other)
        return(gn)
    }else{
      if(modelNum==0){
      	X = NULL
      }else{
      	X=unlist(sapply(effNums,my.design,fixed = TRUE,env = env))
      	X = matrix(X,nrow=N)
      }
      dataRandom = get("dataRandom",env)
      if(!is.null(dataRandom)){
      	
        randomX = apply(data.frame(dataRandom),2,design.mat.single,reduce=FALSE)
        randomX = matrix(randomX,nrow=N)
      }else{
        randomX = NULL
      }
      return(cbind(rep(1,N),X,randomX))
    }
  }
  
}

other.design = function(effNum,fixed=TRUE,env,type="g")
{
  if(effNum==0){
      return(NULL)
  }
  data = get("dataFixed",envir=env)
  #dataRandom = get("dataRandom",envir=env)
  interaction = binary(effNum)$dicotomy
  my.names = colnames(data)[which(interaction)]
  if(sum(interaction)==1){
      df = nlevels(data[,which(interaction)])-fixed
  }else{
      df = prod(unlist(lapply(data[,which(interaction)],nlevels))-fixed)
  }
  if(type=="g"){
    return(df)
  }else if(type=="n"){
    my.names <- gsub(":", "-", my.names, fixed=TRUE)
    return(paste(my.names,collapse=":"))
  }else if(type=='c'){
    if(sum(interaction)==1){
      nlev = nlevels(data[,which(interaction)])
      return(fixedFromRandomProjection(nlev))
    }else{
      eff1 = 2^(which(interaction)[1]-1)
      eff2 = sum(2^(which(interaction)[-1]-1))
      S = other.design(eff1,fixed,env,type="c") %x%  other.design(eff2,fixed,env,type="c")
      return(S)
    }    
  }
}

my.design = function(effNum,fixed=TRUE,env)
{
  if(effNum==0){
    return(NULL)
  }
  X = env$designMatrices[[effNum]]
  if(is.null(X)){
    data = get("dataFixed",envir=env)
    interaction = binary(effNum)$dicotomy
    if(sum(interaction)==1){
      X = design.mat.single(data[,log2(effNum)+1],reduce=fixed)
    }else{
      eff1 = 2^(which(interaction)[1]-1)
      eff2 = sum(2^(which(interaction)[-1]-1))
      X = design.mat.int(my.design(eff1,fixed,env),my.design(eff2,fixed,env))
    }
    env$designMatrices[[effNum]] = X
    return(X)
  }else{
    return(X)
  }
}

# This computes all the model numbers for a given "level"
# (that is, main effects, two-way, three way)
# assuming all effects at lower levels are included
# and no models at higher levels are included
makeModelsSingleLevel <- function(level,nFac)
{
  levels <- choose(nFac,1:nFac)
  totals <- cumsum(levels)
  first <- c()
  last <- c()
  if(level>1){
    first <- rep(TRUE,totals[level-1])    
  }
  if(level<nFac){
    last <- rep(FALSE, totals[nFac] - totals[level])
  }
  #middle <- sapply(0:(2^levels[level]-1),function(n,dim) BayesFactor:::binary(n,dim=dim)$dicotomy,dim=levels[level])
  middle <- as.matrix(expand.grid(rep(list(c(FALSE,TRUE)),levels[level])),ncol=levels[level])
  modelBin <- apply(middle,1,function(v,first,last) c(first,v,last), first=first,last=last)
  apply(modelBin,2,function(v) sum(2^(0:(length(v)-1))[v]))
}

# This computes all models numbers for all levels, as above 
makeModelsAllLevels <- function(nFac){
  modList = unique(unlist(lapply(1:nFac,makeModelsSingleLevel,nFac=nFac)))
  modList = modList[modList!=0]
}

fixedFromRandomProjection <- function(nlevRandom){
  centering=diag(nlevRandom)-(1/nlevRandom)
  S=(eigen(centering)$vectors)[,1:(nlevRandom-1)]
  return(matrix(S,nrow=nlevRandom))
}

