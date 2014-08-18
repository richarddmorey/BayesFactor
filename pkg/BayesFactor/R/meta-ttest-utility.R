meta.ttest.tstat <- function(t,n1,n2=NULL,nullInterval=NULL,rscale)
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
    return(meta.t.bf(t,n,nu,interval=nullInterval,rscale=rscale))
  }
}


meta.t.bf <- function(t,N,df,interval=NULL,rscale){
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
  return(
    list(
      bf = c(null=bf$bf,alt=bf.compl[[1]]),
      properror = c(null=bf$properr,alt=bf.compl[[2]])
    ))
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


meta.t.Metrop <- function(t,n1,n2=NULL,iterations=10000,nullInterval=NULL,rscale,progress=options()$BFprogress, noSample=FALSE, callback = NULL){
  if(length(t)!=length(n1)) stop("lengths of t and n1 must be equal.")
  if(!is.null(n2)){
    if(length(t) != length(n2)) stop("If n2 is defined, it must have the same length as t.")
  }
  iterations = as.integer(iterations)
  
  if(is.null(n2)){
    nu = n1 - 1  
    n <- n1
  }else{
    nu = n1 + n2 - 2
    n <- n1 * n2 / (n1 + n2)
  }
  
  if(progress & !noSample){
    progress = round(iterations/100)
    pb = txtProgressBar(min = 0, max = as.integer(iterations), style = 3) 
  }else{ 
    pb=NULL 
  }
  pbFun = function(samps){ if(progress) setTxtProgressBar(pb, samps)}
  
  
  delta.est = t/sqrt(n)
  mean.delta = sum((delta.est * n)/sum(n))

  if(noSample){
    chains = matrix(NA,1,2)
  }else{
    
    if(is.null(nullInterval)){
      Ubounds = c(0,1)
    }else{
      Ubounds = pnorm(range(nullInterval),mean.delta,sqrt(1/sum(n)))
    }
    
    chains = 1:(iterations+1)
    chains[1] = mean.delta
    
    
    ## Independence Metropolis-Hastings
    for(m in 2:(iterations+1)){
      
      # Cancel the analysis if the callback returns null
      if(is.function(callback))
        if(callback(progress =  m / (iterations + 1) * 1000))
          stop("Operation cancelled.")
      
      candidate = qnorm(runif(1,Ubounds[1],Ubounds[2]),mean.delta,sqrt(1/sum(n)))
    
    ## Metropolis-Hastings acceptance
      z = meta.t.like(candidate,t,n,nu,rscale,log=TRUE) - 
        meta.t.like(chains[m-1],t,n,nu,rscale,log=TRUE) +
        dnorm(chains[m-1],mean.delta,sqrt(1/sum(n)),log=TRUE) - 
        dnorm(candidate,mean.delta,sqrt(1/sum(n)),log=TRUE)
    
      chains[m] = ifelse(rexp(1)>(-z), candidate, chains[m-1])

      if(!(m%%(iterations%/%100))) pbFun(m)
    }
    chains = matrix(chains[-1],iterations)
  }
  
  if(inherits(pb,"txtProgressBar")) close(pb);
  
  acc.rate = mean(diff(chains) != 0)
  message("Independent-candidate M-H acceptance rate: ",round(100*acc.rate),"%")
  return(mcmc(data.frame(delta=chains)))
}



