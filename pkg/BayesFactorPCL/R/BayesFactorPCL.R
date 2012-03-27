logMeanExpLogs = function(v)
{
	N = length(v)
	.Call("RLogMeanExpLogs", as.numeric(v), N, package="BayesFactorPCL")
}

logCumMeanExpLogs = function(v)
{
	N = length(v)
	.Call("RLogCumMeanExpLogs", as.numeric(v), N, package="BayesFactorPCL")
}

nWayAOV.MC = function(y,X,struc,iterations=10000,rscale=1,progress=FALSE,samples=FALSE, gibi=NULL){
	
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
				as.integer(progress), pbFun, new.env(), package="BayesFactorPCL")

	if(progress) close(pb);
	
	#if((returnList[[1]] - nullLike)==Inf){
	#	browser()
	#}
	
	#bf = log(sum(exp(returnList[[2]]-log(iterations)-nullLike)))
	bf = returnList[[1]] - nullLike
	
	if(samples)
	{
		#return(list(returnList[[1]] - nullLike,returnList[[2]]))
		return(list(bf,returnList[[2]]))
	}else{
		#return(returnList[[1]] - nullLike)
		return(bf)
	}
}

ttest.Gibbs = function(y,iterations=10000,rscale=1,null.interval=NULL,progress=TRUE){
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
				progress, pbFun, new.env(), package="BayesFactorPCL")

	if(progress) close(pb);
	priorDens = 1/(pi*rscale)
	postDens = mean(chains[5,])
	BF = postDens/priorDens
	priorArea = pcauchy(interval[2],scale=rscale) - pcauchy(interval[1],scale=rscale)
	postArea = mean(chains[6,])
	BFarea = postArea/priorArea
	
	rownames(chains) = c("mu","sig2","g","delta","CMDE","areaPost")			
	return(list(chains=mcmc(t(chains)),BF=BF,BFarea=BFarea))
}

eqVariance.Gibbs = function(y,iterations=1000,lambda=0.1, M2scale=0.2, sig2.metrop.sd=1 ,tau.metrop.sd=1, M2.metrop.scale=1, newtonSteps=6, progress=TRUE, whichModel=2){
	if(!(whichModel%in%c(1,2))) stop("Argument whichModel model must be either 1 or 2.")
	alpha=0.5
	beta=M2scale^2/2
	N = as.integer(colSums(!is.na(y)))
	J=as.integer(dim(y)[2])
	I=as.integer(dim(y)[1])
	iterations = as.integer(iterations)
	if(progress){
		progress = round(iterations/100)
		pb = txtProgressBar(min = 0, max = as.integer(iterations), style = 3) 
	}else{ 
		pb=NULL 
	}

    pbFun = function(samps){ if(progress) setTxtProgressBar(pb, samps)}
	if(whichModel==1){
		chains = .Call("RgibbsEqVariance", y, N, J, I, lambda, iterations, sig2.metrop.sd, tau.metrop.sd,
				progress, pbFun, new.env(), package="BayesFactorPCL")
		if(progress) close(pb);
		rownames(chains) = c(paste("mu",1:J,sep=""),"CMDE","sig2","sig2.acc","tau","tau.acc")			
		tau.acc = mean(chains[J+5,])
		sig2.acc = mean(chains[J+3,])
		cat("Acceptance rates:\n sig2:",sig2.acc,", tau:",tau.acc,"\n")

		priorDens = integrate(dtau.eqVar,lower=0,upper=Inf,g=rep(1,J),log=FALSE,lambda=lambda)[[1]]
		postDens = mean(chains[J+1,])
		BF = postDens/priorDens
	
		returnList = list(chains=mcmc(t(chains)),BF=BF,acc.rates=c(sig2=sig2.acc,tau=tau.acc))
	}
	if(whichModel==2){
		CI= var(y[,1])*(N[1]-1)/qchisq(c(0.025,0.975),N[1]-1)
		g.metrop.sd = -diff(log(CI))/3 * M2.metrop.scale
		decorr.metrop.sd = g.metrop.sd
		output = .Call("RgibbsEqVarianceM2", y, N, J, I, alpha, beta, iterations, g.metrop.sd, decorr.metrop.sd, as.integer(newtonSteps),
				progress, pbFun, new.env(), package="BayesFactorPCL")
		if(progress) close(pb);
		chains = output[[1]]
		rownames(chains) = c(paste("mu",1:J,sep=""),paste("g",1:J,sep=""),"IWMDE","sig2","sig2g","decorr.acc",paste("g.acc",1:J,sep=""));			
		decorr.acc = mean(chains[2*J+4,])
		g.acc = rowMeans(chains[2*J+ 4 + 1:J,])
		cat("Acceptance rates:\n decorr:",decorr.acc,", g:",g.acc,"\n")

		priorDens = (2*pi)^(-J/2)/gamma(alpha) * gamma(J/2+alpha) * beta^(-J/2)
		postDens = mean(exp(chains[2*J+1,]))
		BF = postDens/priorDens
		returnList = list(chains=mcmc(t(chains)),BF=BF,acc.rates=c(decorr=decorr.acc,g=g.acc))
	}
	
	
	return(returnList)
}

oneWayAOV.Gibbs = function(y,iterations=10000,rscale=1, progress=TRUE, gibi=NULL){
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
				progress, pbFun, new.env(), package="BayesFactorPCL")
	
	if(progress & is.null(gibi)) close(pb);
	rownames(output[[1]]) = c("mu",paste("alpha",1:J,sep=""),"CMDESingle","CMDEDouble","sig2","g")			
	names(output[[2]])=c("logCMDESingle","logCMDEDouble","logCMDESingleKahan","logCMDEDoubleKahan")
	
	logPriorDensDouble = dmvnorm(rep(0,J),rep(0,J),diag(J),log=TRUE)  
	
	#postDensDouble = mean(exp(output[[1]][1+J+2,]))
	logPostDensDouble = logMeanExpLogs(output[[1]][1+J+2,])
	logBFDouble = logPostDensDouble - logPriorDensDouble
	BFDouble = exp(logBFDouble)

	#BFDouble2log = output[[2]][2] - priorDensDouble
	#BFDouble3log = output[[2]][4] - priorDensDouble
	#BFDouble2 = exp(BFDouble2log)
	
	#priorDensSingle = dmvt(rep(0,J),rep(0,J),rscale^2*diag(J),df=1,log=TRUE)
	#postDensSingle = mean(exp(output[[1]][1+J+1,]))
	#BFSingle = postDensSingle/exp(priorDensSingle)
	#BFSingle2log = output[[2]][1] - priorDensSingle
	#BFSingle2 = exp(BFSingle2log)
	#BFSingle3log = output[[2]][3] - priorDensSingle
	
	return(list(chains=mcmc(t(output[[1]])),
			BF=BFDouble,logBF=logBFDouble))
			#bayesFactorRegular=c(single=BFSingle,double=BFDouble),
			#bayesFactorAddLog=c(single=BFSingle2,double=BFDouble2),
			#logBFAddLog=c(single=BFSingle2log,double=BFDouble2log),
			#logBFKahan=c(kahanSingleLogBF=BFSingle3log,kahanDoubleLogBF=BFDouble3log),
			#logCMDE=output[[2]],
			#debug=output[[3]]
			#))
}


marginal.g.oneWay = function(g,F,N,J,rscale)
{
	dfs = (J-1)/(N*J-J)
	omega = (1+(N*g/(dfs*F+1)))/(N*g+1)
	m = log(rscale) - 0.5*log(2*pi) - 1.5*log(g) - rscale^2/(2*g) - (J-1)/2*log(N*g+1) - (N*J-1)/2*log(omega)
	exp(m)
}

oneWayAOV.Quad = function(F,N,J,rscale=1)
{
	1/integrate(marginal.g.oneWay,lower=0,upper=Inf,F=F,N=N,J=J,rscale=rscale)[[1]]
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

ttest.Quad=function(t,n1,n2=0,rscale=1,prior.cauchy=TRUE)
{
	nu=ifelse(n2==0 | is.null(n2),n1-1,n1+n2-2)
	n=ifelse(n2==0 | is.null(n2),n1,(n1*n2)/(n1+n2))
	r2=rscale^2
	marg.like.0=(1+t^2/(nu))^(-(nu+1)/2)
	marg.like.1=ifelse(prior.cauchy,
	integrate(t.joint,lower=0,upper=Inf,t=t,n=n,nu=nu,r2=r2)$value,
		(1+n*r2)^(-.5)*(1+t^2/((1+n*r2)*(nu)))^(-(nu+1)/2))
	return(marg.like.0/marg.like.1)
}


dtau.eqVar.1 = function(x, g, lambda=1, log=FALSE)
{
	if(x<0) return(ifelse(log,-Inf,0))
	J = length(g)
	logd = -J*lgamma(x) - x * (-J*log(x) + sum(g - log(g)) + lambda) + log(lambda)
	#logd = dgamma(g,x,rate=x,log=TRUE) + dexp(x,rate=lambda)
	ifelse(log, logd, exp(logd))
}

dtau.eqVar=Vectorize(dtau.eqVar.1,"x")

integrand.regression=function(g,N,p,R2)
{
       a=.5*((N-p-1)*log(1+g)-(N-1)*log(1+g*(1-R2)))
       exp(a)*dinvgamma(g,shape=.5,scale=N/2)
}

linearReg.Quad=function(N,p,R2) {
	h=integrate(integrand.regression,lower=0,upper=Inf,N=N,p=p,R2=R2)
	return(h$value)
}
