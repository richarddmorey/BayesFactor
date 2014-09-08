
makeTtestHypothesisNames = function(rscale, nullInterval=NULL, mu = 0){
  if(is.null(nullInterval)){
    shortName = paste("Alt., r=",round(rscale,3),sep="")
    longName = paste("Alternative, r = ",rscale,", mu =/= ",mu, sep="")
  }else{
    if(!is.null(attr(nullInterval,"complement"))){
      shortName = paste("Alt., r=",round(rscale,3)," !(",nullInterval[1],"<d<",nullInterval[2], ")",sep="")
      longName = paste("Alternative, r = ",rscale,", mu =/= ",mu, " !(",nullInterval[1],"<d<",nullInterval[2],")",sep="")
    }else{
      shortName = paste("Alt., r=",round(rscale,3)," ",nullInterval[1],"<d<",nullInterval[2],sep="")
      longName = paste("Alternative, r = ",rscale,", mu =/= ",mu, " ",nullInterval[1],"<d<",nullInterval[2],sep="")
    }
  }
  return(list(shortName=shortName,longName=longName))
}


ttestBF_oneSample = function(x, mu, nullInterval, rscale, posterior, callback, ... ){
  rscale = rpriorValues("ttestOne",,rscale)
  hypNames = makeTtestHypothesisNames(rscale, nullInterval, mu = mu)
  
  mod1 = BFoneSample(type = "JZS", 
                        identifier = list(formula = "y ~ 1", nullInterval = nullInterval), 
                        prior=list(rscale=rscale, mu=mu, nullInterval = nullInterval),
                        shortName = hypNames$shortName,
                        longName = hypNames$longName)
  if(posterior)
    return(posterior(mod1, data = data.frame(y=x), callback = callback, ...))
  
  bf1 = compare(numerator = mod1, data = data.frame(y=x))
  
  if(!is.null(nullInterval)){
    mod2 = mod1
    attr(mod2@identifier$nullInterval, "complement") = TRUE
    attr(mod2@prior$nullInterval, "complement") = TRUE
    hypNames = makeTtestHypothesisNames(rscale, mod2@identifier$nullInterval, mu = mu)
    mod2@shortName = hypNames$shortName
    mod2@longName = hypNames$longName
    
    bf2 = compare(numerator = mod2, data = data.frame(y=x))    
    return(c(bf1,bf2))
  }else{
    return(bf1)
  }    
  
}

ttestBF_indepSample = function(formula, data, mu, nullInterval, rscale, posterior, callback, ... ){
  checkFormula(formula, data, analysis = "indept")
  
  rscale = rpriorValues("ttestTwo",,rscale)
  hypNames = makeTtestHypothesisNames(rscale, nullInterval, mu = mu)
  
  mod1 = BFindepSample(type = "JZS", 
                          identifier = list(formula = stringFromFormula(formula), nullInterval = nullInterval), 
                          prior=list(rscale=rscale, mu=mu, nullInterval = nullInterval),
                          shortName = hypNames$shortName,
                          longName = hypNames$longName
  )
  
  if(posterior)
    return(posterior(mod1, data = data, callback = callback, ...))
  
  bf1 = compare(numerator = mod1, data = data)
  
  if(!is.null(nullInterval)){
    mod2 = mod1
    attr(mod2@identifier$nullInterval, "complement") = TRUE
    attr(mod2@prior$nullInterval, "complement") = TRUE
    hypNames = makeTtestHypothesisNames(rscale, mod2@identifier$nullInterval, mu = mu)
    mod2@shortName = hypNames$shortName
    mod2@longName = hypNames$longName
    
    bf2 = compare(numerator = mod2, data = data)
    return(c(bf1, bf2))
  }else{
    return(bf1)
  }  
  
  
}  


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
  if(is.null(callback)) callback=function(...) as.integer(0)

  if(noSample){
    chains = matrix(NA,6,2)
  }else{
    chains = .Call("RgibbsOneSample", as.numeric(y), n, rscale, iterations, do.interval, as.numeric(interval),
                 progress, pbFun, callback, new.env(), package="BayesFactor")
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


