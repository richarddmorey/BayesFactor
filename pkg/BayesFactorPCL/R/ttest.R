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
	lbfarea = log(postArea) - log(1-postArea) - (log(priorArea) - log(1-priorArea))
	
	rownames(chains) = c("mu","sig2","g","delta","CMDE","areaPost")			
	if(logbf){
		return(list(chains=mcmc(t(chains)),BF=-lbf,BFarea=-lbfarea))
	}else{
		return(list(chains=mcmc(t(chains)),BF=exp(-lbf),BFarea=exp(-lbfarea)))
	}
}

t.joint=function(g,t,n,nu,r2)
{
	t1=-.5*log(1+n*g*r2)
	t2=(-(nu+1)/2)*log(1+t^2/((1+n*g*r2)*(nu)))
	return(dinvgamma(g,.5,.5)*exp(t1+t2))
}

