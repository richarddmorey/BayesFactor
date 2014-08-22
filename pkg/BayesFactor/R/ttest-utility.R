ttest.Gibbs = function(y=NULL,t=NULL,n=NULL,iterations=10000,rscale="medium",nullInterval=NULL,progress=options()$BFprogress,logbf=FALSE,noSample=FALSE, callback=NULL){
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
  
  rscale = rpriorValues("ttestOne",,rscale)
  iterations = as.integer(iterations)
  if(progress & !noSample){
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
  
  if(noSample){
    chains = matrix(NA,6,2)
  }else{
    chains = .Call("RgibbsOneSample", as.numeric(y), n, rscale, iterations, do.interval, as.numeric(interval),
                 progress, pbFun, new.env(), package="BayesFactor")
  }
  
  if(inherits(pb,"txtProgressBar")) close(pb);
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


