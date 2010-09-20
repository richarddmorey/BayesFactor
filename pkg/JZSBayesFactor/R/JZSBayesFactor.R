JZSgibbs = function(y,iterations=1000,rscale=1,progress=TRUE){
	N = as.integer(length(y))
	iterations = as.integer(iterations)
	if(progress){
		progress = round(iterations/100)
		pb = txtProgressBar(min = 0, max = as.integer(iterations), style = 3) 
	}else{ 
		pb=NULL 
	}

    pbFun = function(samps){ if(progress) setTxtProgressBar(pb, samps)}
	
	chains = .Call("RgibbsOneSample", y, N, rscale, iterations,
				progress, pbFun, new.env(), package="JZSBayesFactor")

	if(progress) close(pb);
	priorDens = 1/(pi*rscale)
	postDens = mean(chains[5,])
	BF = postDens/priorDens
	rownames(chains) = c("mu","sig2","g","delta","CMDE")			
	return(list(chains=mcmc(t(chains)),bayesFactor=BF))
}

eqVarGibbs = function(y,iterations=1000,lambda=1, sig2.metrop.sd=1 ,tau.metrop.sd=1, progress=TRUE){
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

	chains = .Call("RgibbsEqVariance", y, N, J, I, lambda, iterations, sig2.metrop.sd, tau.metrop.sd,
				progress, pbFun, new.env(), package="JZSBayesFactor")

	if(progress) close(pb);
	rownames(chains) = c(paste("mu",1:J,sep=""),"CMDE","sig2","sig2.acc","tau","tau.acc")			
	tau.acc = mean(chains[J+5,])
	sig2.acc = mean(chains[J+3,])
	cat("Acceptance rates:\n sig2:",sig2.acc,", tau:",tau.acc,"\n")

	priorDens = integrate(dtau.eqVar,lower=0,upper=Inf,g=rep(1,J),log=FALSE,lambda=lambda)[[1]]
	postDens = mean(chains[J+1,])
	BF = postDens/priorDens
	
	return(list(chains=mcmc(t(chains)),bayesFactor=BF))

}

oneWayAOVGibbs = function(y,iterations=1000,rscale=1, progress=TRUE){
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

	chains = .Call("RgibbsOneWayAnova", y, N, J, I, rscale, iterations,
				progress, pbFun, new.env(), package="JZSBayesFactor")

	if(progress) close(pb);
	rownames(chains) = c("mu",paste("beta",1:J,sep=""),"CMDESingle","CMDEDouble","sig2","g")			

	priorDensDouble = dmvnorm(rep(0,J),rep(0,J),diag(J))  
	postDensDouble = mean(chains[1+J+2,])
	BFDouble = postDensDouble/priorDensDouble
	
	priorDensSingle = dmvt(rep(0,J),rep(0,J),rscale^2*diag(J),df=1,log=FALSE)
	postDensSingle = mean(chains[1+J+1,])
	BFSingle = postDensSingle/priorDensSingle
	
	return(list(chains=mcmc(t(chains)),bayesFactor=c(single=BFSingle,double=BFDouble)))

}


marginal.g.oneWay = function(g,F,N,J,rscale)
{
dfs = (J-1)/(N*J-J)
omega = (1+(N*g/(dfs*F+1)))/(N*g+1)
m = log(rscale) - 0.5*log(2*pi) - 1.5*log(g) - rscale^2/(2*g) - (J-1)/2*log(N*g+1) - (N*J-1)/2*log(omega)
exp(m)
}

F.value.bf = function(F,N,J,rscale=1)
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

t.value.bf=function(t,n1,n2=0,sd=1,prior.cauchy=T)
{
nu=ifelse(n2==0 | is.null(n2),n1-1,n1+n2-2)
n=ifelse(n2==0 | is.null(n2),n1,(n1*n2)/(n1+n2))
r2=sd*sd
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

