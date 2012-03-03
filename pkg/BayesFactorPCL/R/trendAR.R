trendtest.Gibbs.AR = function(before,after,iterations=1000,r.scaleInt=1,r.scaleSlp=1,alphaTheta=1,betaTheta=5,sdMet=.3, progress=TRUE,return.chains=FALSE)
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
	
	chains = .Call("RgibbsTwoSampleAR_trend", y, N, X, ncol(X), r.scaleInt, r.scaleSlp, alphaTheta, betaTheta, iterations, sdMet, progress, pbFun, new.env(), package="BayesFactorPCL")
	
	if(progress) close(pb)
	
	#dim(out[[2]]) = c(,iterations)
	dim(chains) = c(ncol(X) + 7, iterations)
	chains = data.frame(t(chains))
	colnames(chains) = c(paste("beta",1:ncol(X),sep=""),"sig2","g1","g2","rho","densFullRes", "densSlpRes","densIntRes")
	
	ldens = apply(chains[,9:11],2,logMeanExpLogs)
	nulllogdens = c(
				dcauchy(0,log=TRUE) - log(r.scaleSlp) + 
				dcauchy(0,log=TRUE) - log(r.scaleInt),
				dcauchy(0,log=TRUE) - log(r.scaleSlp),
				dcauchy(0,log=TRUE) - log(r.scaleInt)
				)
	logbf = ldens - nulllogdens
	
	acc = mean(diff(chains$rho)!=0)
	cat("\n","rho acceptance rate:",acc,"\n")
	
	if(return.chains)
	{
		return(list(logbf=logbf,chains=mcmc(chains),acc=acc,debug=NULL))
	}else{
		return(c(logbf=logbf))
	}
}
