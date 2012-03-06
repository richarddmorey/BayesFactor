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
