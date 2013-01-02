
oneWayAOV.Gibbs = function(y,iterations=10000,rscale="medium", progress=TRUE, gibi=NULL, logbf=FALSE){
	
  rscale = rpriorValues("allNways","fixed",rscale)
  N = as.integer(colSums(!is.na(y)))
	J=as.integer(dim(y)[2])
	I=as.integer(dim(y)[1])
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
	
	output = .Call("RgibbsOneWayAnova", y, N, J, I, rscale, iterations,
				progress, pbFun, new.env(), package="BayesFactor")
	
	if(progress & is.null(gibi)) close(pb);
	rownames(output[[1]]) = c("mu",paste("beta",1:J,sep=""),"CMDESingle","CMDEDouble","sig2","g")			
	names(output[[2]])=c("logCMDESingle","logCMDEDouble","logCMDESingleKahan","logCMDEDoubleKahan")
	
	logPriorDensDouble = dmvnorm(rep(0,J),rep(0,J),diag(J),log=TRUE)  
	
	logPostDensDouble = logMeanExpLogs(output[[1]][1+J+2,])
	lbf = logPostDensDouble - logPriorDensDouble

	if(logbf){
		return(list(chains=mcmc(t(output[[1]])), BF=-lbf))
	}else{
		return(list(chains=mcmc(t(output[[1]])), BF=exp(-lbf)))
	}
	

}

oneWayAOV.Quad = function(F,N,J,rscale="medium",logbf=FALSE,error.est=FALSE)
{
  rscale = rpriorValues("allNways","fixed",rscale)
  integral = integrate(marginal.g.oneWay,lower=0,upper=Inf,F=F,N=N,J=J,rscale=rscale)
  properror = exp(log(integral[[2]]) - log(integral[[1]]))
	bf = ifelse(logbf, log(integral[[1]]), integral[[1]])
  if(error.est){
		return(c(bf=bf, properror=properror))
	}else{
		return(bf)
	}
}


marginal.g.oneWay = function(g,F,N,J,rscale)
{
	dfs = (J-1)/(N*J-J)
	omega = (1+(N*g/(dfs*F+1)))/(N*g+1)
	m = log(rscale) - 0.5*log(2*pi) - 1.5*log(g) - rscale^2/(2*g) - (J-1)/2*log(N*g+1) - (N*J-1)/2*log(omega)
	exp(m)
}