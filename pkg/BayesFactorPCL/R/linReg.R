linearReg.Quad=function(N,p,R2,rscale=1,logbf=FALSE) {
	h=integrate(integrand.regression,lower=0,upper=Inf,N=N,p=p,R2=R2,rscaleSqr=rscale^2)
	if(logbf){
		return(log(h$value))
	}else{
		return(h$value)
	}
}


integrand.regression=function(g,N,p,R2,rscaleSqr=1)
{
       a=.5*((N-p-1)*log(1+g)-(N-1)*log(1+g*(1-R2)))
       exp(a)*dinvgamma(g,shape=.5,scale=rscaleSqr*N/2)
}

linearReg.Gibbs <- function(y, covariates, iterations = 10000, rscale = 1, progress = TRUE, gibi=NULL){
  	X <- apply(covariates,2,function(v) v - mean(v))
  	y = matrix(y,ncol=1)
  	N = length(y)
  	p = ncol(X)
  	Cn = diag(N) - matrix(1,N,N)/N 
  	XtX = t(X)%*%X
  	XtCnX = t(X)%*%Cn%*%X
  	XtCny = t(X)%*%Cn%*%y
  	Cny = Cn%*%y

    sig2start = sum( (X%*%solve(XtCnX)%*%XtCny - Cny)^2 ) / N
      
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

	chains = .Call("RGibbsLinearReg", 
						 as.integer(iterations), 
						 Cny,
						 X,
						 XtX,
						 XtCnX,
						 XtCny,
						 as.integer(N),
						 as.integer(p),
						 rscale,
             sig2start,
						 progress,
						 pbFun,
						 new.env(),
						 package="BayesFactor")
						 
  if(progress & is.null(gibi)) close(pb);
  chains = t(chains)
  
  colnames(chains) = c(paste("beta",1:p,sep=""),"sig2","g")
  return(mcmc(chains))

}
