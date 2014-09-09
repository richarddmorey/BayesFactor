makeMetaTtestHypothesisNames = function(rscale, nullInterval=NULL){
  if(is.null(nullInterval)){
    shortName = paste("Alt., r=",round(rscale,3),sep="")
    longName = paste("Alternative, r = ",rscale,", delta =/= 0", sep="")
  }else{
    if(!is.null(attr(nullInterval,"complement"))){
      shortName = paste("Alt., r=",round(rscale,3)," (",nullInterval[1],"<d<",nullInterval[2],")",sep="")
      longName = paste("Alternative, r = ",rscale,", delta =/= 0 !(",nullInterval[1],"<d<",nullInterval[2],")",sep="")
    }else{
      shortName = paste("Alt., r=",round(rscale,3)," ",nullInterval[1],"<d<",nullInterval[2],sep="")
      longName = paste("Alternative, r = ",rscale,", delta =/= 0 ",nullInterval[1],"<d<",nullInterval[2],sep="")
    }
  }
  return(list(shortName=shortName,longName=longName))
}


meta.ttest.tstat <- function(t,n1,n2=NULL,nullInterval=NULL,rscale, complement = FALSE)
{
  
  if( ((length(n1) != length(n2)) & !is.null(n2)) | (length(n1) != length(t))){
    stop("Number of t statistics must equal number of sample sizes.")  
  }
  if(!is.null(n2)){
    rscale = rpriorValues("ttestTwo",,rscale)
  }else{
    rscale = rpriorValues("ttestOne",,rscale)
  }
  if(is.null(n2)){
    nu = n1 - 1  
    n <- n1
  }else{
    nu = n1 + n2 - 2
    n <- n1 * n2 / (n1 + n2)
  }
  
  if(is.null(nullInterval)){
    return(meta.t.bf(t,n,nu,rscale=rscale))
  }else{
    return(meta.t.bf(t,n,nu,interval=nullInterval,rscale=rscale, complement = complement))
  }
}


meta.t.bf <- function(t,N,df,interval=NULL,rscale, complement = FALSE){
  if(length(interval)!=2 & !is.null(interval))
    stop("argument interval must have two elements.")
  
  if(is.null(interval)){
    return(meta.bf.interval(-Inf,Inf,t,N,df,rscale))
  }
  
  interval = unique(sort(interval))
  
  if(interval[1]==-Inf & interval[2]==Inf){
    return(meta.bf.interval(-Inf,Inf,t,N,df,rscale))
  }
  
  if(any(is.infinite(interval))){
    bf = meta.bf.interval(interval[1],interval[2],t,N,df,rscale)
    if(interval[1]==-Inf){
      bf.compl = meta.bf.interval(interval[2],Inf,t,N,df,rscale)  
    }else{
      bf.compl = meta.bf.interval(-Inf,interval[1],t,N,df,rscale)
    } 
  }else{
    logPriorProbs = pcauchy(c(-Inf,interval,Inf),scale=rscale,log.p=TRUE)
    prior.interval1 = logExpAminusExpB(logPriorProbs[2], logPriorProbs[1])
    prior.interval3 = logExpAminusExpB(logPriorProbs[4], logPriorProbs[3])
    
    prior.interval.1.3 = logMeanExpLogs(c(prior.interval1,prior.interval3)) + log(2)
    
    bf1 = meta.bf.interval(-Inf,interval[1],t,N,df,rscale)
    bf = meta.bf.interval(interval[1],interval[2],t,N,df,rscale)
    bf3 = meta.bf.interval(interval[2],Inf,t,N,df,rscale)
    
    bf.compl = sumWithPropErr(bf1$bf + prior.interval1,
                              bf3$bf + prior.interval3,
                              bf1$properror,
                              bf3$properror)
    bf.compl[1] = bf.compl[1] - prior.interval.1.3
  }
  
  if(complement){
    return(
      list(
        bf = bf.compl[[1]],
        properror = bf.compl[[2]]
      ))
  }else{
    return(
      list(
        bf = bf$bf,
        properror = bf$properr
      ))
  }
}


meta.bf.interval <- function(lower,upper,t,N,df,rscale=1){
  nullLike = sum(dt(t,df,log=TRUE))
  logPriorProbs = pcauchy(c(upper,lower),scale=rscale,log.p=TRUE)
  prior.interval = logExpAminusExpB(logPriorProbs[1], logPriorProbs[2])
  delta.est = t/sqrt(N)
  mean.delta = sum((delta.est * N)/sum(N))
  log.const = meta.t.like(mean.delta,t,N,df,rscale,log=TRUE)
  intgl = integrate(meta.t.like,lower,upper,t=t,N=N,df=df,rscale=rscale,log.const=log.const)
  
  val = log(intgl[[1]]) + log.const - prior.interval - nullLike
  err = exp(log(intgl[[2]]) - val)
  
  return(
    list(
      bf = val,
      properror = err
    )
  )
}

meta.t.like <- Vectorize(function(delta,t,N,df,rscale=1,log.const=0,log=FALSE){
  ans = suppressWarnings(
    sum(dt(t,df,ncp=delta*sqrt(N),log=TRUE)) + dcauchy(delta,scale=rscale,log=TRUE) - log.const
  )
  if(log){
    return(ans)
  }else{
    return(exp(ans))
  }
},"delta")


meta.t.Metrop <- function(t, n1, n2=NULL, nullModel, iterations=10000, nullInterval=NULL, rscale, progress=options()$BFprogress, noSample=FALSE, callback = NULL, callbackInterval = 1){
  if(length(t)!=length(n1)) stop("lengths of t and n1 must be equal.")
  if(!is.null(n2)){
    if(length(t) != length(n2)) stop("If n2 is defined, it must have the same length as t.")
  }
  iterations = as.integer(iterations)
  
  if( is.null(n2) ){
    n2 = n1*0
    twoSample = FALSE
  }else{
    twoSample = TRUE
  }
  
  progress = as.logical(progress)
  if(is.null(callback) | !is.function(callback)) callback=function(...) as.integer(0)
  
  if(is.null(nullInterval) | nullModel){
    doInterval = FALSE
    nullInterval = c(-Inf, Inf)
    intervalCompl = FALSE
  }else{
    doInterval = TRUE
    intervalCompl = ifelse(!is.null(attr(nullInterval,"complement")),TRUE,FALSE)
    nullInterval = range(nullInterval)
  }
  
  if(noSample){
    chains = matrix(as.numeric(NA),1,1)
  }else{
    if(nullModel) rscale = 0
    chains = metropMetaTRcpp(t, n1, n2, twoSample, rscale, iterations, 
                           doInterval, nullInterval, intervalCompl, nullModel, 
                           progress, callback, callbackInterval)
    if(!nullModel){
      acc.rate = mean(diff(chains) != 0)
      message("Independent-candidate M-H acceptance rate: ",round(100*acc.rate),"%")
    }
  }
  
  return(mcmc(data.frame(delta=chains)))
}



