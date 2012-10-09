logMeanExpLogs = function(v)
{
	N = length(v)
	.Call("RLogMeanExpLogs", as.numeric(v), N, package="BayesFactor")
}

logCumMeanExpLogs = function(v)
{
	N = length(v)
	.Call("RLogCumMeanExpLogs", as.numeric(v), N, package="BayesFactor")
}

testNwaymLike = function(g,y,Xm)
{
	P = dim(Xm)[2]
	N = dim(Xm)[1]
	
	y = matrix(y,N)
	Xc = t(t(Xm) - colMeans(Xm))
	
	m0 = -0.5*(N - 1) * log(var(y)*(N-1))
	
	.Call("RjeffmlikeNWayAov",
		t(Xc)%*%Xc,
		t(Xc)%*%y,
		var(y)*(N-1),
		N, P, g, 
		package="BayesFactor") - m0
}

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

ttest.Quad=function(t,n1,n2=0,rscale=sqrt(2)/2,prior.cauchy=TRUE,logbf=FALSE)
{
	nu=ifelse(n2==0 | is.null(n2),n1-1,n1+n2-2)
	n=ifelse(n2==0 | is.null(n2),n1,(n1*n2)/(n1+n2))
	r2=rscale^2
	marg.like.0=(1+t^2/(nu))^(-(nu+1)/2)
	marg.like.1=ifelse(prior.cauchy,
	integrate(t.joint,lower=0,upper=Inf,t=t,n=n,nu=nu,r2=r2)$value,
		(1+n*r2)^(-.5)*(1+t^2/((1+n*r2)*(nu)))^(-(nu+1)/2))
	lbf = log(marg.like.0) - log(marg.like.1)
	if(logbf){
		return(-lbf)
	}else{
		return(exp(-lbf))
	}
}

ttest.Gibbs = function(y,iterations=10000,rscale=sqrt(2)/2,null.interval=NULL,progress=TRUE,logbf=FALSE){
	N = as.integer(length(y))
	iterations = as.integer(iterations)
	if(progress){
		progress = round(iterations/100)
		pb = txtProgressBar(min = 0, max = as.integer(iterations), style = 3) 
	}else{ 
		pb=NULL 
	}
	
	if(is.null(null.interval)){
		do.interval=0
		interval = c(-Inf,Inf)
	}else{
		if(length(null.interval)!=2){
			stop("null.interval must be a vector of length 2.")
		}
		do.interval=1
		interval=sort(as.numeric(null.interval))
	}
    
    pbFun = function(samps){ if(progress) setTxtProgressBar(pb, samps)}
	
	chains = .Call("RgibbsOneSample", y, N, rscale, iterations, do.interval, interval,
				progress, pbFun, new.env(), package="BayesFactor")

	if(progress) close(pb);
	priorDens = 1/(pi*rscale)
	postDens = mean(chains[5,])
	lbf = log(postDens) - log(priorDens)
	priorArea = pcauchy(interval[2],scale=rscale) - pcauchy(interval[1],scale=rscale)
	postArea = mean(chains[6,])
	lbfarea = log(postArea) - log(priorArea)
	
	rownames(chains) = c("mu","sig2","g","delta","CMDE","areaPost")			
	if(logbf){
		return(list(chains=mcmc(t(chains)),BF=-lbf,BFarea=-lbfarea))
	}else{
		return(list(chains=mcmc(t(chains)),BF=exp(-lbf),BFarea=exp(-lbfarea)))
	}
}


oneWayAOV.Gibbs = function(y,iterations=10000,rscale=1/2, progress=TRUE, gibi=NULL, logbf=FALSE){
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
	rownames(output[[1]]) = c("mu",paste("alpha",1:J,sep=""),"CMDESingle","CMDEDouble","sig2","g")			
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

oneWayAOV.Quad = function(F,N,J,rscale=1/2,logbf=FALSE)
{
	lbf = -log(integrate(marginal.g.oneWay,lower=0,upper=Inf,F=F,N=N,J=J,rscale=rscale)[[1]])
	if(logbf){
		return(-lbf)
	}else{
		return(exp(-lbf))
	}
}


marginal.g.oneWay = function(g,F,N,J,rscale)
{
	dfs = (J-1)/(N*J-J)
	omega = (1+(N*g/(dfs*F+1)))/(N*g+1)
	m = log(rscale) - 0.5*log(2*pi) - 1.5*log(g) - rscale^2/(2*g) - (J-1)/2*log(N*g+1) - (N*J-1)/2*log(omega)
	exp(m)
}


dinvgamma = function (x, shape, scale = 1) 
{
    if (shape <= 0 | scale <= 0) {
        stop("Shape or scale parameter negative in dinvgamma().\n")
    }
    alpha <- shape
    beta <- scale
    log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 
        1) * log(x) - (beta/x)
    return(exp(log.density))
}

t.joint=function(g,t,n,nu,r2)
{
	t1=-.5*log(1+n*g*r2)
	t2=(-(nu+1)/2)*log(1+t^2/((1+n*g*r2)*(nu)))
	return(dinvgamma(g,.5,.5)*exp(t1+t2))
}

integrand.regression=function(g,N,p,R2)
{
       a=.5*((N-p-1)*log(1+g)-(N-1)*log(1+g*(1-R2)))
       exp(a)*dinvgamma(g,shape=.5,scale=N/2)
}

linearReg.Quad=function(N,p,R2,logbf=FALSE) {
	h=integrate(integrand.regression,lower=0,upper=Inf,N=N,p=p,R2=R2)
	if(logbf){
		return(log(h$value))
	}else{
		return(h$value)
	}
}
