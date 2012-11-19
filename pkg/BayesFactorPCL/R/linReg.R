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
