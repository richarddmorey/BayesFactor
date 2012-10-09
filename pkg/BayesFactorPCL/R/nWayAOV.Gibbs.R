
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

