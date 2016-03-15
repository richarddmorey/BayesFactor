makePropHypothesisNames = function(rscale, nullInterval=NULL, p){
  if(is.null(nullInterval)){
    shortName = paste("Alt., p0=",p,", r=",round(rscale,3),sep="")
    longName = paste("Alternative, p0 = ",p,", r = ",rscale,", p =/= p0", sep="")
  }else{
    if(!is.null(attr(nullInterval,"complement"))){
      shortName = paste("Alt., p0=",p,", r=",round(rscale,3)," !(",nullInterval[1],"<p<",nullInterval[2],")",sep="")
      longName = paste("Alternative, p0 = ", p ,", r = ",rscale,", p =/= p0 !(",nullInterval[1],"<p<",nullInterval[2],")",sep="")
    }else{
      shortName = paste("Alt., p0=",p,", r=",round(rscale,3)," ",nullInterval[1],"<p<",nullInterval[2],sep="")
      longName = paste("Alternative, p0 = ",p,", r = ",rscale,", p =/= p0 ",nullInterval[1],"<p<",nullInterval[2],sep="")
    }
  }
  return(list(shortName=shortName,longName=longName))
}


like.prop.test = Vectorize(function(lo, y, N, p, rscale, log = FALSE, log.const = 0, shift = 0, scale = 1){
  ans = sum(dbinom(y, N, plogis(lo*scale + shift), log=TRUE)) + dlogis(lo*scale + shift, qlogis(p), rscale, log=TRUE) - log.const
  ifelse(log, ans, exp(ans))
}
,"lo")


prop.test.bf.interval <- function(y, N, p, rscale, nullInterval){ 
  intervalProb = plogis(nullInterval, qlogis(p), rscale, log.p = TRUE)
  prior.interval = logExpXminusExpY(intervalProb[2], intervalProb[1])
  
  beta.a = sum(y) + 1
  beta.b = sum(N) - sum(y) + 1
  lo.est = qlogis(beta.a/(beta.a + beta.b))
  p.var = beta.a * beta.b / ( ( beta.a + beta.b )^2 * ( beta.a + beta.b + 1 ) )
  lo.sd = sqrt( p.var / dlogis( lo.est )^2 )
  
  log.const = like.prop.test(lo.est, y = y, N = N, p = p, rscale = rscale, log = TRUE)
  
  intgl = integrate(like.prop.test, (nullInterval[1] - lo.est)/lo.sd, (nullInterval[2] - lo.est)/lo.sd, y = y, N = N, p = p, rscale = rscale, log.const = log.const, shift = lo.est, scale = lo.sd)
  nullLike = sum(dbinom(y, N, p, log = TRUE))
  
  val = log(intgl[[1]]*lo.sd) + log.const - prior.interval - nullLike
  err = exp(log(intgl[[2]]) + log(lo.sd) - val)
  
  return(
    list(
      bf = val,
      properror = err,
      method = "quadrature"
    )
  )
}

prop.test.bf <- function(y, N, p, rscale, interval, complement){
  if(length(interval)!=2 & !is.null(interval))
    stop("argument interval must have two elements.")
  
  if(is.null(interval)){
    return(prop.test.bf.interval(y, N, p, rscale, c(-Inf,Inf)))
  }

  interval = range(interval)
  interval = qlogis(interval)
  
  if(interval[1]==-Inf & interval[2]==Inf){
    if(complement){
      return(list(bf=NA,properror=NA)) 
    }else{
      return(prop.test.bf.interval(y, N, p, rscale, interval))
    }
  }
  
  if(any(is.infinite(interval))){
    bf = prop.test.bf.interval(y, N, p, rscale, interval)
    if(interval[1]==-Inf){
      bf.compl = prop.test.bf.interval(y, N, p, rscale, c(interval[2], Inf))
    }else{
      bf.compl = prop.test.bf.interval(y, N, p, rscale, c(-Inf,interval[1]))
    } 
  }else{
    logPriorProbs = plogis(c(-Inf,interval,Inf), qlogis(p), scale=rscale,log.p=TRUE)
    prior.interval1 = logExpXminusExpY(logPriorProbs[2], logPriorProbs[1])
    prior.interval3 = logExpXminusExpY(logPriorProbs[4], logPriorProbs[3])
    
    prior.interval.1.3 = logMeanExpLogs(c(prior.interval1,prior.interval3)) + log(2)
    
    bf1 = prop.test.bf.interval(y, N, p, rscale, c(-Inf,interval[1]))
    bf = prop.test.bf.interval(y, N, p, rscale, interval)
    bf3 = prop.test.bf.interval(y, N, p, rscale, c(interval[2],Inf))
    
    bf.compl = sumWithPropErr(bf1[['bf']] + prior.interval1,
                              bf3[['bf']] + prior.interval3,
                              bf1[['properror']],
                              bf3[['properror']])
    bf.compl[1] = bf.compl[1] - prior.interval.1.3
  }
  
  if(complement){
    return(
      list(
        bf = bf.compl[[1]],
        properror = bf.compl[[2]],
        method = "quadrature"
      ))
  }else{
    return(
      list(
        bf = bf[['bf']],
        properror = bf[['properror']],
        method = bf[['method']]
      ))
  }
}



proportion.Metrop <- function(y, N, nullModel, iterations=10000, nullInterval=NULL, p, rscale, progress=options()$BFprogress, noSample=FALSE, callback = NULL, callbackInterval = 1){
  if(length(y)!=length(N)) stop("lengths of t and n1 must be equal.")
  iterations = as.integer(iterations)
   
  progress = as.logical(progress)
  if(is.null(callback) | !is.function(callback)) callback=function(...) as.integer(0)
  
  if(is.null(nullInterval) | nullModel){
    doInterval = FALSE
    nullInterval = c(0, 1)
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
    chains = metropProportionRcpp(y, N, p, rscale, iterations, doInterval, 
                         qlogis(nullInterval), intervalCompl, nullModel, 
                         progress, callback, callbackInterval) 
    if(!nullModel & !noSample){
      acc.rate = mean(diff(chains) != 0)
      message("Independent-candidate M-H acceptance rate: ",round(100*acc.rate),"%")
    }
  }
  chains = mcmc(data.frame(logodds = chains, p = plogis(chains)))
  return(chains)
}


