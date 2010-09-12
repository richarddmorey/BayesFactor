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
	
	rownames(chains) = c("mu","sig2","g","delta","CMDE")			
	return(mcmc(t(chains)))
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
	return(mcmc(t(chains)))
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
