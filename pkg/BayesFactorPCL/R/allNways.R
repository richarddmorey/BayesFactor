allNways = function(y,dataFixed=NULL,dataRandom=NULL,iterations = 10000, which.models="withmain", progress=TRUE, rscaleFixed=.5, rscaleRandom=1, logbf=FALSE, multicore=FALSE, ...)
{
  if( !(which.models %in% c("all","top","withmain"))){
    stop("Invalid value for which.models: must be 'all', 'top', or 'withmain'")
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
  
  nFac = dim(dataFixed)[2]
  if(nFac==1) which.models='all'
  bfEnv = new.env(parent = baseenv())
  designs = list()
  designs[[2^nFac]] = matrix(nrow=0,ncol=0)

  bfEnv$designMatrices = designs
  bfEnv$dataFixed = dataFixed
  bfEnv$y = y
  bfEnv$totalN = length(as.vector(y))
  bfEnv$dataRandom = dataRandom

  
  if(multicore){
    allResults <- all.Nways.env.mc(env=bfEnv,iterations=iterations, which.models, progress=FALSE, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom, ...)
  }else{
    allResults <- all.Nways.env(env=bfEnv,iterations=iterations, which.models, progress=progress, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom, ...)
  }
  bfs = as.numeric(allResults[1,])
  names(bfs) = allResults[2,]
  bfs = c(null=0,bfs)
  	
  if(!is.null(dataRandom))
  {
  	nullMod = as.numeric(
                nWayAOV2(0,bfEnv,iterations=iterations, which.models, 
                      rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom, ...)[1]
                )
  	bfs = bfs - nullMod
	bfs[1] = 0
  }
  if(which.models=="top"){
  	topModel = ((2^(2^nFac-1))-1)
  	topModelName= paste(joined.design(topModel,env=bfEnv,other="n"),collapse=" + ")
  	topModelIndex = which(names(bfs)==topModelName)
  	bfs = bfs - bfs[topModelIndex]
  }
  if(logbf){
  	return(sort(bfs, decreasing=TRUE))
  }else{
  	return(sort(exp(bfs), decreasing=TRUE))
  }
}

all.Nways.env = function(env, which.models, progress=progress, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom,...){
  data = env$dataFixed
	nFac = dim(data)[2]
	topMod = ((2^(2^nFac-1))-1)
	if(which.models=="all"){
		modNums = 1:topMod
	}else if(which.models=="top"){
		nDig = 2^nFac-1
		mods <- c(colSums((1-diag(nDig))*2^(0:(nDig-1))),topMod)
		modNums <- as.list(mods)
	}else if( which.models=="withmain"){
    modNums = makeModelsAllLevels(nFac)
	}
  if(progress){
    results <- pbsapply(modNums,nWayAOV2,env=env, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom,...)     
  }else{
    results <- sapply(modNums,nWayAOV2,env=env, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom,...)
  }
	return(results)
}


#### Multi core version
all.Nways.env.mc = function(env, which.models,progress=progress, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom,...){
  if(!require(doMC)){
    warning("Required package (doMC) missing for multicore functionality. Falling back to single core functionality.")
    allResults <- all.Nways.env(env=env, which.models=which.models, progress=progress, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom, ...)
    return(allResults)
  }  
  
  registerDoMC()
  if(getDoParWorkers()==1){
    warning("Multicore specified, but only using 1 core. Set options(cores) to something >1.")
  }
  
  data = env$dataFixed
  nFac = dim(data)[2]
  topMod = ((2^(2^nFac-1))-1)
  if(which.models=="all"){
    modNums = 1:topMod
  }else if(which.models=="top"){
    nDig = 2^nFac-1
    mods <- c(colSums((1-diag(nDig))*2^(0:(nDig-1))),topMod)
    modNums <- as.list(mods)
  }else if( which.models=="withmain"){
    modNums = makeModelsAllLevels(nFac)
  }
  
  # Taken from http://stackoverflow.com/questions/10984556/is-there-way-to-track-progress-on-a-mclapply
  # Multicore progress bar: does not work!
  if(progress){
    bfs <- local({
      f <- fifo(tempfile(), open="w+b", blocking=TRUE)
      if (inherits(fork(), "masterProcess")) {
        # Child
        progressSoFar <- 0.0
        cat(progressSoFar,"\n")
        while (progressSoFar < 1 & !isIncomplete(f)) {
          msg <- readBin(f, "double")
          progressSoFar <- progressSoFar + as.numeric(msg)
          cat(sprintf("Progress: %.2f%%\n", progressSoFar * 100))
        } 
        exit()
      }
      numJobs <- length(modNums)
      progressCallback = function(){
        writeBin(1/numJobs, f)
      }
    
      bfs <- foreach(gIndex=modNums,.combine='cbind', .options.multicore=mcoptions) %dopar% nWayAOV2(gIndex,env=env, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom, progressCallback = progressCallback,...) 
    
      close(f)
      bfs
    })
  }else{
    # No progress bar
    bfs <- foreach(gIndex=modNums,.combine='cbind', .options.multicore=mcoptions) %dopar% nWayAOV2(gIndex,env=env, rscaleFixed=rscaleFixed, rscaleRandom=rscaleRandom,...) 
  }  
  return(bfs)
}

buildModelInfo <- function(modNum, env) {
  X = joined.design(modNum,env=env)
  my.name = paste(joined.design(modNum,env=env,other="n"),collapse=" + ")
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
    names=my.name
    ))
}


nWayAOV2 = function(modNum,env, rscaleFixed, rscaleRandom, progressCallback=NULL, ...)
{
  modInfo = buildModelInfo(modNum, env)
  y = env$y
  
  bfs = c(nWayAOV.MC(y, modInfo$X, modInfo$struc, samples=FALSE,logbf=TRUE, progress=FALSE, 
  						rscale = c(rscaleFixed + modInfo$g.groups*0, rscaleRandom=rscaleRandom + modInfo$gr.groups*0),...), 
  						modInfo$names)
  names(bfs)=modInfo$names
  
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

do.row.mult = function(v,p1,p2)
{
  v1 = v[1:p1]
  v2 = v[p1 + 1:p2]
  as.vector(outer(v1,v2,'*'))
}

joined.design = function(modelNum, env, other=NULL)
{
  N = get("totalN",env)
  dataRandom = get("dataRandom",env)
  if(modelNum==0 & is.null(dataRandom)){
    #if(!is.null(other)){
    #  if(other=="g"){
    #  	return(NULL)
    #  }else if(other=="n"){
    #  	return("<null>")
    #  }
    #}else{
      return(matrix(1,nrow=N))
    #}
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

