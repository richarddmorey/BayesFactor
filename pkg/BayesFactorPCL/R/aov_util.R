# Taken from the WLE package source by Claudio Agostinelli <claudio at unive.it>

binary <- function(x, dim) {

   if (x==0) {
       pos <- 1
   } else {
       pos <- floor(log(x, 2))+1
   }

   if (!missing(dim)) {
       if (pos<=dim) {
           pos <- dim
       } else {
           warning("the value of `dim` is too small")
       }  
   }

   bin <- rep(0, pos)
   dicotomy <- rep(FALSE, pos)
   for (i in pos:1) {
        bin[i] <- floor(x/2^(i-1))
        dicotomy[i] <- bin[i]==1
        x <- x-((2^(i-1))*bin[i])
   }
   return(list(binary=bin, dicotomy=dicotomy))
}



design.mat.single=function(v,reduce=TRUE)
{
  nlev = nlevels(as.factor(v))
  v = as.integer(as.factor(v))
  X = matrix(0,nrow=length(v),ncol=nlev)
  X[1:length(v) + length(v)*(v-1)] = 1
  if(reduce){
    centering=diag(nlev)-(1/nlev)
    S=(eigen(centering)$vectors)[,1:(nlev-1)]
    return(X%*%S)
  }else{
    return(X)
  }
}

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
  dataRandom = get("dataRandom",envir=env)
  interaction = binary(effNum)$dicotomy
  my.names = colnames(data)[which(interaction)]
  if(sum(interaction)==1){
      df = nlevels(data[,which(interaction)])-1
  }else{
      df = prod(unlist(lapply(data[,which(interaction)],nlevels))-1)
  }
  if(type=="g"){
    return(df)
  }else if(type=="n"){
    return(paste(my.names,collapse=":"))
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

nWayAOV2 = function(modNum,env,samples=FALSE, logFunction = cat,...)
{
  #logFunction(paste(modNum,"\n"))
  flush.console()
  X = joined.design(modNum,env=env)
  y = env$y
  g.groups = joined.design(modNum,env=env,other="g")
  my.name = paste(joined.design(modNum,env=env,other="n"),collapse=" + ")
  dataRandom = env$dataRandom
    if(!is.null(dataRandom)){
    gr.groups = unlist(lapply(data.frame(dataRandom),nlevels))
  }else{
    gr.groups = NULL
  } 
  bfs = c(nWayAOV.MC(y,X,c(g.groups,gr.groups),samples=samples,...),my.name)
  if(!samples)
  { 
    names(bfs)=my.name
  }
  return(bfs)
}

allNways = function(y,dataFixed=NULL,dataRandom=NULL,iterations = 1000, samples=FALSE, only.top=FALSE,...)
{
  nFac = dim(dataFixed)[2]
  if(nFac==1) only.top=FALSE
  bfEnv = new.env(parent = baseenv())
  designs = list()
  designs[[2^nFac]] = matrix(nrow=0,ncol=0)

  bfEnv$designMatrices = designs
  bfEnv$dataFixed = dataFixed
  bfEnv$y = y
  bfEnv$totalN = length(as.vector(y))
  bfEnv$dataRandom = dataRandom
  if(!samples){
  	allResults <- all.Nways.env(env=bfEnv,iterations=iterations, samples=FALSE, only.top,...)
	bfs = as.numeric(allResults[1,])
	names(bfs) = allResults[2,]
  	bfs = c(null=0,bfs)
  	
  	if(!is.null(dataRandom))
  	{
  		nullMod = as.numeric(nWayAOV2(0,bfEnv,iterations=iterations, samples=FALSE, only.top,...)[1])
  		bfs = bfs - nullMod
		bfs[1] = 0
  	}
  	return(sort(bfs))
  }else{
  	allResults <- all.Nways.env(env=bfEnv,iterations=iterations, samples=TRUE, only.top,...)
  	bfs <- unlist(lapply(allResults,function(lst) lst[[1]] ))
  	bfs = c(null=0,bfs)
  	if(!is.null(dataRandom))
  	{
  		nullMod = nWayAOV2(0,bfEnv,iterations=iterations, samples=TRUE,...)
  		bfs = (bfs - nullMod[[1]])
		bfs[1] = 0
		nullSamp = nullMod[[2]]
  	}else{
  		nullSamp = 0
  	}
 	
 	my.names <- unlist(lapply(allResults,function(lst) lst[[3]] ))
  	my.names = c("null",my.names)
	names(bfs)=my.names

	samples = c(null=list(nullSamp),lapply(allResults,function(lst) lst[[2]] )) 

  	return(list(bfs,samples))
  }

}

all.Nways.env = function(env,samples, only.top=FALSE,...){
  data = env$dataFixed
	nFac = dim(data)[2]
	topMod = ((2^(2^nFac-1))-1)
	if(!only.top){
		modNums = 1:topMod
	}else{
		nDig = 2^nFac-1
		mods <- c(colSums((1-diag(nDig))*2^(0:(nDig-1))),topMod)
		modNums <- as.list(mods)
	}
	if(!samples){
		return(sapply(modNums,nWayAOV2,env=env,samples=FALSE,...))
	}else{
		return(lapply(modNums,nWayAOV2,env=env,samples=TRUE,...))
		
	}
}

