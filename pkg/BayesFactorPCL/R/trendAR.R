trendtest.Gibbs.AR = function(before, after, iterations=1000, intArea=c(-.2,.2), slpArea=c(-.2,.2), r.scaleInt=1, r.scaleSlp=1,alphaTheta=1,betaTheta=5,sdMet=.3, progress=TRUE,return.chains=FALSE)
{


	y = c(before,after)
	N = length(y)	
	
	iterations = as.integer(iterations)

	if(progress){
		progress = round(iterations/100)
		pb = txtProgressBar(min = 0, max = as.integer(iterations), style = 3) 
	}else{ 
		pb=NULL 
	}
	
    pbFun = function(samps){ if(progress) setTxtProgressBar(pb, samps)}

	treat = c(rep(-0.5,length(before)),rep(0.5,length(after)))
	timePoint = 1:N - length(before) - .5
	X = cbind(1,treat)
	X = cbind(X,X*timePoint)
	
	out = .Call("RgibbsTwoSampleAR_trend", y, N, X, ncol(X), r.scaleInt, r.scaleSlp, alphaTheta, betaTheta, intArea[1], intArea[2], slpArea[1], slpArea[2], iterations, sdMet, progress, pbFun, new.env(), package="BayesFactorPCL")
	
	if(progress) close(pb)
	
	#dim(out[[2]]) = c(,iterations)
	dim(out[[1]]) = c(ncol(X) + 4, iterations)
	out[[1]] = data.frame(t(out[[1]]))
	colnames(out[[1]]) = c(paste("beta",1:ncol(X),sep=""),"sig2","g1","g2","rho")#,"densFullRes", "densSlpRes","densIntRes","areaInt","areaSlp")
	
	ldens = out[[2]][1:3] - log(iterations)#apply(chains[,9:11],2,logMeanExpLogs)
	nulllogdens = c(
				dcauchy(0,log=TRUE) - log(r.scaleSlp) + 
				dcauchy(0,log=TRUE) - log(r.scaleInt),
				dcauchy(0,log=TRUE) - log(r.scaleSlp),
				dcauchy(0,log=TRUE) - log(r.scaleInt)
				)
	logbf = ldens - nulllogdens
	
	areas = log(out[[2]][4:5]) - log(iterations)#log(colMeans(chains[,12:13]))
	nullAreas = log(c(
				diff(pcauchy(intArea,scale=r.scaleInt)),
				diff(pcauchy(slpArea,scale=r.scaleSlp))
				))
	acc = mean(diff(out[[1]]$rho)!=0)
	cat("\n","rho acceptance rate:",acc,"\n")
	
	if(return.chains)
	{
		return(list(logbf=logbf, chains=mcmc(out[[1]]), acc=acc, logbfArea=areas - nullAreas,debug=NULL))
	}else{
		return(c(logbf=logbf))
	}
}

trendtest.MC.AR = function(before, after, iterations=1000, r.scaleInt=1, r.scaleSlp=1,alphaTheta=1,betaTheta=5, progress=TRUE,return.chains=FALSE)
{
	y = c(before,after)
	N = length(y)	
	
	iterations = as.integer(iterations)

	if(progress){
		progress = round(iterations/100)
		pb = txtProgressBar(min = 0, max = as.integer(iterations), style = 3) 
	}else{ 
		pb=NULL 
	}
	
    pbFun = function(samps){ if(progress) setTxtProgressBar(pb, samps)}

	treat = c(rep(-0.5,length(before)),rep(0.5,length(after)))
	timePoint = 1:N - length(before) - .5
	X0 = cbind(1,treat)
	X1 = X0*timePoint
	X = cbind(X0,X1)
	
	nullLike = log(trendtest.nullMargLikeAR(y,X[,c(1,3)],alphaTheta,betaTheta))

	out = .Call("MCAR_trend", y, N, alphaTheta, betaTheta, r.scaleInt^2, r.scaleSlp^2, X[,c(1,3)], X[,c(2,4)], iterations, progress, pbFun, new.env(), package="BayesFactorPCL")

	if(progress) close(pb)
	
	altLike = out[[1]]
	
	bfs = c(nullLike-altLike[1],altLike[2]-altLike[1],altLike[3]-altLike[1])
	
	if(return.chains){	
		return(list(bfs,nullLike,altLike,t(out[[2]])))
	}else{
		return(list(bfs,nullLike,altLike))
	}


}	

trendtest.nullMargLikeAR = function(y,X0,alphaTheta=1,betaTheta=5){
	N = length(y)
	fun = Vectorize(trendtest.nullMargLikeAR.theta, "theta")
	integrate(fun,0,1,y=y,X0=X0,alphaTheta=alphaTheta,betaTheta=betaTheta,N=N)[[1]]
}

trendtest.nullMargLikeAR.theta = function(theta,y,N=length(y),X0,alphaTheta=1,betaTheta=5){
	ret = .Call("MCnullMargLogLikeAR_trend",theta,y,N,alphaTheta,betaTheta,X0,package="BayesFactorPCL")
	return(exp(ret))
}

trendtest.altMargLikeAR.thetag1g2 = function(theta,g1,g2,y,N=length(y),X0,X1,rInt=1, rSlp=1,alphaTheta=1,betaTheta=5){
	ret = .Call("MCmargLogLikeAR_trendR",theta,g1,g2,y,N,alphaTheta,betaTheta,rInt,rSlp,X0,X1,package="BayesFactorPCL")
	return(ret)
}