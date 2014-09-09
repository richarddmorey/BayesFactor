
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


ttestOneSample.Gibbs = function(y, nullModel, iterations, rscale, nullInterval, progress=options()$BFprogress, noSample=FALSE, callback=NULL, callbackInterval = 1){
  n = as.integer(length(y))
  
  rscale = ifelse(nullModel,1,rpriorValues("ttestOne",,rscale))
  iterations = as.integer(iterations)
  progress = as.logical(progress)

  if(is.null(nullInterval)){
    do.interval = FALSE
    nullInterval = c(-Inf,Inf)
    complement = FALSE
  }else{
    if(length(nullInterval)!=2){
      stop("nullInterval must be a vector of length 2.")
    }
    do.interval=TRUE
    complement = ifelse(!is.null(attr(nullInterval,"complement")),TRUE,FALSE)
  }
  
  if(is.null(callback) | !is.function(callback)) callback=function(...) as.integer(0)

  if(noSample){
    chains = matrix(as.numeric(NA),4,2)
  }else{
    chains = gibbsOneSampleRcpp(mean(y), var(y), n, rscale, iterations, do.interval, nullInterval, complement, nullModel, progress, callback, callbackInterval) 
  }
    
  colnames(chains) = c("mu","sig2","delta","g")			
  if(nullModel){
    return(mcmc(chains[,-4]))
  }else{
    return(mcmc(chains))    
  }
}

ttestIndepSample.Gibbs = function(formula, data, nullModel, iterations, rscale, nullInterval, progress=options()$BFprogress, noSample=FALSE, callback=NULL, callbackInterval = 1){
  
  depVar = as.character(formula[[2]])
  
  if(nullModel){
    # If sampling from the null model, it doesn't matter how we divide the data
    fakeIndepVar = c(1,1,rep(2, nrow(data) - 2))
    ybar = tapply(data[,depVar], fakeIndepVar, mean)   
    s2 = tapply(data[,depVar], fakeIndepVar, var)
    N =  tapply(data[,depVar], fakeIndepVar, length)
    rscale = 1
  }else{
    indepVar = as.character(formula[[3]])
    spltData = rev(split(data[,depVar], factor(data[,indepVar])))
    ybar = sapply(spltData, mean)
    
    if(length(ybar)!=2) stop("Incorrect number of levels in independent variable.")
    
    s2 = sapply(spltData, var)
    N =  sapply(spltData, length)
    grp.names = names(spltData)
    effect.direction = paste(rev(grp.names), collapse = " - ")
    rscale = rpriorValues("ttestTwo",,rscale)
  }
  
  
  iterations = as.integer(iterations)
  progress = as.logical(progress)
  
  if(is.null(nullInterval)){
    do.interval = FALSE
    nullInterval = c(-Inf,Inf)
    complement = FALSE
  }else{
    if(length(nullInterval)!=2){
      stop("nullInterval must be a vector of length 2.")
    }
    do.interval=TRUE
    complement = ifelse(!is.null(attr(nullInterval,"complement")),TRUE,FALSE)
  }
  
  if(is.null(callback) | !is.function(callback)) callback=function(...) as.integer(0)
  
  if(noSample){
    chains = matrix(as.numeric(NA),5,2)
  }else{
    chains = gibbsTwoSampleRcpp(ybar, s2, N, rscale, iterations, do.interval, nullInterval, complement, nullModel, progress, callback, callbackInterval)
  }
  
  
  if(nullModel){
    colnames(chains) = c("mu","beta","sig2","delta","g")    	
    return(mcmc(chains[,-5]))
  }else{
    colnames(chains) = c("mu",paste("beta (",effect.direction,")",sep=""),"sig2","delta","g")    	
    return(mcmc(chains))    
  }
}



