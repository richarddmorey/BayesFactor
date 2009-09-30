#Jeff Rouder
#Last Edit 9/08

#Richard Morey
#Edited April 4/09


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

ols.g=function(g,Rsqr,n,v)
{
t1=((n-2)/2)*log(1+g*v)
t2=(-(n-1)/2)*log(1+(1-Rsqr)*g*v)
return(dinvgamma(g,.5,n/2)*exp(t1+t2))
}

ols.value.bf=function(F,n,v=1,prior.cauchy=T)
{
Rsqr=F/(F+n-2)
bfinv=ifelse(prior.cauchy,
integrate(ols.g,lower=0,upper=Inf,Rsqr=Rsqr,n=n,v=v)$value,
exp(((n-2)/2)*log(1+n*v)+(-(n-1)/2)*log(1+(1-Rsqr)*n*v)))
return(1/bfinv)
}

#############################


PCLBF=function(t,N1,N2=0,prior.cauchy=TRUE,null.d=.2,pi0=.5,scale=1,MCMC=FALSE,CIconf=.95,...)
{
BF=new("PCLBF")

BF@t.value=t
BF@N1=N1
BF@N2=N2
BF@prior.cauchy=prior.cauchy
BF@MCMC=MCMC
BF@pi0=pi0
BF@null.d=null.d
BF@prior.scale=scale

t.df=ifelse(N2==0 | is.null(N2),N1-1,N1+N2-2)
my.N=ifelse(N2==0 | is.null(N2),N1,(N1*N2)/(N1+N2))
BF@num.samples=ifelse(N2==0 | is.null(N2),1,2)


nu=ifelse(prior.cauchy,1,Inf)
pi2=pt(abs(null.d)/scale,nu)-pt(-abs(null.d)/scale,nu)

BF@JZS=t.value.bf(t,N1,N2,prior.cauchy=prior.cauchy,sd=scale)
area=area.bf(t,N1,N2,nu,null.d,MCMC=MCMC,r.scale=scale,...)


BF@NOH=area$bf.area
BF@marginal.d=area$post.d
BF@MCMC.iters=area$iterations
BF@MCMC.95CI=area$bfMCMC95
BF@MCMC.acceptance=c(area$acc.rate.delta,area$acc.rate.eta2)
BF@MCMC.chains=area$chains
BF@freq.CI=freqCIdelta(t,N1,alpha=1-CIconf)
BF@freq.conf=CIconf
BF@freq.p = 2*(1-pt(abs(t),t.df))
BF@hybrid = pi0*pi2*BF@JZS*BF@NOH + pi0*BF@JZS - pi0*pi2*BF@JZS + BF@NOH - pi0*BF@NOH



return(BF)
}


hybrid.bf=function(t,N,pi0=.5,prior.cauchy=TRUE,null.d=.2,MCMC=FALSE,scale=1)
{
nu=ifelse(prior.cauchy,1,Inf)
pi2=pt(abs(null.d)/scale,nu)-pt(-abs(null.d)/scale,nu)
B1 = t.value.bf(t,N,prior.cauchy=prior.cauchy)
B2 = area.bf(t,N,df.prior=nu,null.d=null.d,MCMC=MCMC)
B3 = pi0*pi2*B1*B2 + pi0*B1 - pi0*pi2*B1 + B2 - pi0*B2
B3
}

area.bf.integrate = function(t, N, df.prior=1, null.d=0.2, MCMC=FALSE,eta2.lo=0,eta2.up=Inf,d.lo=-Inf,d.up=Inf,r.scale=1)
{
	null.d=abs(null.d)

	priorarea = pt(null.d/r.scale,df.prior) - pt(-null.d/r.scale,df.prior)
	priorodds = priorarea/(1-priorarea)

	postarea.null  = integrate(vect.margd.post, -null.d, null.d, t=t, N=N, df.prior=df.prior,r.scale=r.scale,rel.tol=max(50*.Machine$double.eps,0.5e-28))$value
	postarea.total = integrate(vect.margd.post, d.lo, d.up, t=t, N=N, df.prior=df.prior,r.scale=r.scale,rel.tol=max(50*.Machine$double.eps,0.5e-28))$value

	post.d=function(x) vect.margd.post(x, t, N, df.prior=1,eta2.lo=eta2.lo,eta2.up=eta2.up,r.scale=r.scale)/postarea.total


	BF = (postarea.null/(postarea.total - postarea.null))/priorodds
      
      output=list(bf.area=BF,bfMCMC95=0,post.d=post.d,t=t,N=N,nu=df.prior,iterations=0,acc.rate.delta=0,acc.rate.eta2=0,chains=0)

	return(output)
}

area.bf = function(t, N1,N2=NULL, df.prior=1, null.d=0.2, MCMC=FALSE,r.scale=1,...)
{
	if(MCMC){
		out=area.bf.mcmc(t,N1,df.prior,null.d,,r.scale=r.scale,...)
	}else{
		out=area.bf.integrate(t,N1,df.prior,null.d,r.scale=r.scale,...)		
	}
out
}


area.bf.mcmc=function(t,N,nu=1,null.d=.2,iters=10000,metscale=2,pnull=.5,chains=FALSE,r.scale=1){
        null.d=abs(null.d)

        pnull.cor=ifelse(nu==Inf,pnorm(null.d/r.scale)-pnorm(-null.d/r.scale),pt(null.d/r.scale,nu)-pt(-null.d/r.scale,nu))

        area.cor=pnull.cor/(1-pnull.cor)
        prior.cor=pnull/(1-pnull)

        a=N/2
        b=(N-1)*1/2
        sdfullcond.eta2=b/((a-1)*sqrt(a-2))
        sdfullcond.delta=1/sqrt(N)

        metscale.eta2=sdfullcond.eta2*metscale
        metscale.delta=sdfullcond.delta*metscale

        eta2.samp=1:iters*NA
        delta.samp=1:iters*NA

        eta2.samp[1]=1
        delta.samp[1]=t/sqrt(N)


        for(i in 2:iters){
                eta2.samp[i]=samp.eta2.t(eta2.samp[i-1],delta.samp[i-1],t,N,metscale.eta2)
                delta.samp[i]=samp.delta.t(delta.samp[i-1],eta2.samp[i],t,N,metscale.delta,nu=nu,r=r.scale)
        }

        postmean.delta=mean(delta.samp)


        acc.rate.delta=mean(diff(delta.samp)!=0)
        acc.rate.eta2=mean(diff(eta2.samp)!=0)

        postprobnull=mean(delta.samp>-null.d & delta.samp<null.d)
        postoddsnull=postprobnull/(1-postprobnull)

        bf.area=(postoddsnull/area.cor)/prior.cor

        innull=delta.samp>-null.d & delta.samp<null.d


        phat=mean(innull)
        mcmc95=phat+c(-1,1)*1.96*sqrt(phat*(1-phat))/sqrt(effectiveSize(as.numeric(innull)))
        bfCI95=((mcmc95/(1-mcmc95))/area.cor)/prior.cor

        
	  post.d=approxfun(density(delta.samp)$x,density(delta.samp)$y)

        if(chains){
			my.chains=as.mcmc(data.frame(delta.samp,eta2.samp,as.numeric(innull)))
		}else{
			my.chains=0
		}        

        output=list(bf.area=bf.area,bfMCMC95=bfCI95,post.d=post.d,t=t,N=N,nu=nu,iterations=iters,acc.rate.delta=acc.rate.delta,acc.rate.eta2=acc.rate.eta2,chains=my.chains)

  
        output
}



