ttestBF <- function(x, y = NULL, formula = NULL, mu = 0, nullInterval = NULL, 
                    paired = FALSE, data = NULL, rscale="medium"){
  
  rscale = rpriorValues("ttest",,rscale)
  
  if( (is.null(formula) & is.null(y)) | (!is.null(y) & paired) ){
    if(paired){
      # check that the two vectors have same length
      if(length(x)!=length(y)) stop("Length of x and y must be the same if paired=TRUE.")
      x = x - y
    }
    modFull = BFoneSample(type = "JZS one sample", 
                                 identifier = list(formula = "y ~ 1"), 
                                 prior=list(rscale=rscale, mu=mu),
                                 shortName = paste("Alt., r=",round(rscale,3),sep=""),
                                 longName = paste("Alternative, r = ",rscale,", mu =/= ",mu, sep="")
                         )
    if(is.null(nullInterval)){
      bf = compare(numerator = modFull, data = data.frame(y=x))
      return(bf)
    }else{
      nullInterval = range(nullInterval)
      modInterval = BFoneSample(type = "JZS one sample", 
                                 identifier = list(formula = "y ~ 1",nullInterval = nullInterval), 
                                 prior=list(rscale=rscale, mu=mu, nullInterval = nullInterval),
                                 shortName = paste("Alt., r=",round(rscale,3)," ",nullInterval[1],"<d<",nullInterval[2],sep=""),
                                 longName = paste("Alternative, r = ",rscale,", mu =/= ",mu, " ",nullInterval[1],"<d<",nullInterval[2],sep="")
      )      
      bf = compare(numerator = modInterval, data = data.frame(y=x))
      return(bf)
    }
  }else if(!is.null(y) & !paired){
    data = data.frame(y = c(x,y), 
                      group = factor(c(rep("x",length(x)),rep("y",length(y))))
                      )
    formula = y ~ group
  }
  if(!is.null(formula)){ # formula
    if(paired) stop("Cannot use 'paired' with formula.")
    if(is.null(data)) stop("'data' needed for formula.")
    
    ivs = attr(terms(formula),"term.labels")
    if(length(ivs) > 1) stop("Only one independent variable allowed for t test.")
    dataTypes = "fixed"
    names(dataTypes) = ivs
    
    modFull = BFlinearModel(type = "JZS independent samples", 
                          identifier = list(formula = deparse(formula)), 
                          prior=list(rscale=rscale, mu=mu),
                          dataTypes = dataTypes,
                          shortName = paste("Alt., r=",round(rscale,3),sep=""),
                          longName = paste("Alternative, r = ",rscale,", mu1-mu2 =/= ",mu, sep="")
    )
    if(is.null(nullInterval)){
      bf = compare(numerator = modFull, data = data)
      return(bf)
    }else{
      nullInterval = range(nullInterval)
      modInterval = BFlinearModel(type = "JZS independent samples", 
                                identifier = list(formula = deparse(formula),nullInterval = nullInterval), 
                                prior=list(rscale=rscale, mu=mu, nullInterval = nullInterval),
                                dataTypes = dataTypes,
                                shortName = paste("Alt., r=",round(rscale,3)," ",nullInterval[1],"<d<",nullInterval[2],sep=""),
                                longName = paste("Alternative, r = ",rscale,", mu =/= ",mu, " ",nullInterval[1],"<d<",nullInterval[2],sep="")
      )
     bf = compare(numerator = modInterval, data = data)
      return(bf)
    }
    
  }
}

ttest.Quad=function(t,n1,n2=0,nullInterval=NULL,rscale="medium",logbf=FALSE,error.est=FALSE)
{
  rscale = rpriorValues("ttest",,rscale)
  
  nu=ifelse(n2==0 | is.null(n2),n1-1,n1+n2-2)
	n=ifelse(n2==0 | is.null(n2),n1,(n1*n2)/(n1+n2))
	r2=rscale^2
	marg.like.0=(1+t^2/(nu))^(-(nu+1)/2)
  if(is.null(nullInterval)){
    integral = integrate(t.joint,lower=0,upper=Inf,t=t,n=n,nu=nu,r2=r2)  
	  marg.like.1 = integral$value
	  prop.error = integral$abs.error / marg.like.1
    lbf = log(marg.like.1) - log(marg.like.0)
  }else{
    areabf = ttestAreaNull(t, n1, n2, nullInterval=nullInterval, rscale=rscale)
    lbf = areabf$bf
    prop.error = areabf$properror
  }
  if(logbf){
    if(error.est){
      return(list(bf = lbf, properror=prop.error))
    }else{
      return(lbf)
    }  	
  }else{
    if(error.est){
      return(list(bf = exp(lbf), properror=prop.error))
    }else{
      return(exp(lbf))
    }  	
	}
}

ttest.Gibbs = function(y=NULL,t=NULL,n=NULL,iterations=10000,rscale="medium",null.interval=NULL,progress=TRUE,logbf=FALSE){
	if( (is.null(t) | is.null(n)) & !is.null(y) ){
    n = as.integer(length(y))
  }else if(!is.null(t) & !is.null(n)){
    # Create some fake data with needed parameters to pass
    y = rnorm(n)
    y = (y - mean(y))/sd(y)
    y = y + t / sqrt(n)
    n = as.integer(n)
	}else{
    stop("Insufficient data: either t, or both t and n, must be given.")
	}
  rscale = rpriorValues("ttest",,rscale)
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
	
	chains = .Call("RgibbsOneSample", y, n, rscale, iterations, do.interval, interval,
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

ttestAreaNull <- function(t, n1, n2=0, nullInterval=c(-.2,.2), rscale=1, safeInt = .9999)
{
  nu=ifelse(n2==0 | is.null(n2),n1-1,n1+n2-2)
  n=ifelse(n2==0 | is.null(n2),n1,(n1*n2)/(n1+n2))
  
  nullInterval = range(nullInterval)
  safeRange = t/sqrt(n) + c(-1,1) * qt(1-(1-safeInt)/2,nu)/sqrt(n)
  
  priorOdds = diff(pcauchy(nullInterval,scale=rscale))

  nullInterval[1] = max(nullInterval[1],safeRange[1]) 
  nullInterval[2] = min(nullInterval[2],safeRange[2]) 
  
  ifelse(nullInterval[1]<safeRange[1],safeRange[1],nullInterval[1])
  
  allIntegral = integrate(function(delta,tstat,n1,nu,rscale){
    exp(dt(t,df = nu, ncp = delta * sqrt(n1), log=TRUE) + dcauchy(delta, scale=rscale, log=TRUE))
  }, safeRange[1], safeRange[2], tstat=t,n1=n,nu=nu, rscale=rscale)[[1]]
  
  areaIntegral = integrate(function(delta,tstat,n1,nu,rscale,const=1){
    exp(dt(t,df = nu, ncp = delta * sqrt(n1), log=TRUE) + dcauchy(delta, scale=rscale, log=TRUE) - log(const))
  }, nullInterval[1], nullInterval[2], tstat=t,n1=n,nu=nu,rscale=rscale,const=allIntegral)

  
  # encompassing vs point null
  vsNull = ttest.Quad(t, n1, n2, rscale=rscale, logbf=TRUE, error.est=TRUE)
  
  val = areaIntegral[[1]]
  err = areaIntegral[[2]]
  
  err = err / val 
  err = sqrt(err^2 + vsNull[['properror']]^2)
  val = log(val) -  log(priorOdds) + vsNull[['bf']]
  
  complArea = 1-areaIntegral[[1]]
  errCompl = areaIntegral[[2]] / complArea
  errCompl = sqrt(errCompl^2 + vsNull[['properror']]^2)  
  valCompl = log(complArea) - log(1-priorOdds) + vsNull[['bf']]
  
  return(
    list(
      bf = c(null=val,alt=valCompl),
      properror = c(null=err,alt=errCompl)
    ))
}