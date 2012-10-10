nWayAOV.MC = function(y,X,struc,iterations=10000,rscale=1,progress=FALSE,samples=FALSE, gsamples=FALSE, gibi=NULL, logbf=FALSE){
	
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
				as.integer(progress), pbFun, new.env(), package="BayesFactor")

	if(progress) close(pb);
	
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


nWayAOV.Gibbs <- function(y,X,struc,iterations=10000,rscale=1,progress=TRUE)
{
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
	labels = c("mu",paste("beta",1:p,sep="_"),"sig2",paste("g",1:length(struc),sep="_"))
	colnames(chains) = labels
	return(chains)
}

