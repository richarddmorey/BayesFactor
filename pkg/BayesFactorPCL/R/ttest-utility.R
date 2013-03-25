ttest.Gibbs = function(y=NULL,t=NULL,n=NULL,iterations=10000,rscale="medium",nullInterval=NULL,progress=options()$BFprogress,logbf=FALSE){
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
  
  if(is.null(nullInterval)){
    do.interval=0
    interval = c(-Inf,Inf)
  }else{
    if(length(nullInterval)!=2){
      stop("nullInterval must be a vector of length 2.")
    }
    do.interval=1
    interval=sort(as.numeric(nullInterval))
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
  vsNull = ttest.tstat(t, n1, n2, rscale=rscale)
  
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

