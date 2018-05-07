makeCorrHypothesisNames = function(rscale, nullInterval=NULL){
  if(is.null(nullInterval)){
    shortName = paste("Alt., r=",round(rscale,3),sep="")
    longName = paste("Alternative, r = ",rscale,", rho =/= 0", sep="")
  }else{
    if(!is.null(attr(nullInterval,"complement"))){
      shortName = paste("Alt., r=",round(rscale,3)," !(",nullInterval[1],"<rho<",nullInterval[2],")",sep="")
      longName = paste("Alternative, r = ",rscale,", rho =/= 0 !(",nullInterval[1],"<rho<",nullInterval[2],")",sep="")
    }else{
      shortName = paste("Alt., r=",round(rscale,3)," ",nullInterval[1],"<rho<",nullInterval[2],sep="")
      longName = paste("Alternative, r = ",rscale,", rho =/= 0 ",nullInterval[1],"<rho<",nullInterval[2],sep="")
    }
  }
  return(list(shortName=shortName,longName=longName))
}


corr.test.bf.interval <- function(y, x, rscale, nullInterval){
  intervalProb = pbeta(nullInterval/2+.5, 1/rscale, 1/rscale, log.p=TRUE)
  prior.interval = logExpXminusExpY(intervalProb[2], intervalProb[1])

  r = cor(y, x, use="pairwise.complete.obs")
  n = length(x) - sum(is.na(y) | is.na(x))
  # if(n != length(x)) note(paste("Ignored",sum(is.na(y) | is.na(x)),
  #                               "rows containing missing observations."))

  intgl = .bfCorNumerical(n, r, rscale, nullInterval[1], nullInterval[2])
  val = log(intgl$value) - prior.interval
  err = NA

  return(
    list(
      bf = val,
      properror = err,
      method = "mixed"
    )
  )
}

corr.test.bf <- function(y, x, rscale, nullInterval, complement){
  if(length(nullInterval)!=2 & !is.null(nullInterval))
    stop("argument interval must have two elements.")

  r = cor(y, x, use="pairwise.complete.obs")
  n = length(x) - sum(is.na(y) | is.na(x))
  if(n != length(x)) message(paste("Ignored",sum(is.na(y) | is.na(x)),
                                "rows containing missing observations."))
  if(is.null(nullInterval)){
    return(.bf10Exact(n, r, rscale))
  }

  interval = range(nullInterval)

  if(nullInterval[1]<=-1 & nullInterval[2]>=1){
    if(complement){
      return(list(bf=NA,properror=NA))
    }else{
      return(return(.bf10Exact(n, r, rscale)))
    }
  }

  if(any(abs(nullInterval)>=1)){
    bf = corr.test.bf.interval(y, x, rscale, nullInterval)
    if(interval[1]<=-1){
      bf.compl = corr.test.bf.interval(y, x, rscale, c(nullInterval[2], 1))
    }else{
      bf.compl = corr.test.bf.interval(y, x, rscale, c(-1,nullInterval[1]))
    }
  }else{
    logPriorProbs = pbeta(c(-1,nullInterval,1)/2+.5, 1/rscale, 1/rscale, log.p=TRUE)
    prior.interval1 = logExpXminusExpY(logPriorProbs[2], logPriorProbs[1])
    prior.interval3 = logExpXminusExpY(logPriorProbs[4], logPriorProbs[3])

    prior.interval.1.3 = logMeanExpLogs(c(prior.interval1,prior.interval3)) + log(2)

    bf1 = corr.test.bf.interval(y, x, rscale, c(-1,nullInterval[1]))
    bf = corr.test.bf.interval(y, x, rscale, nullInterval)
    bf3 = corr.test.bf.interval(y, x, rscale, c(nullInterval[2],1))

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



correlation.Metrop <- function(y, x, nullModel, iterations=10000, nullInterval=NULL, rscale, progress=getOption('BFprogress', interactive()), noSample=FALSE, callback = NULL, callbackInterval = 1){
  if(length(y)!=length(x)) stop("lengths of y and x must be equal.")
  iterations = as.integer(iterations)

  progress = as.logical(progress)
  if(is.null(callback) | !is.function(callback)) callback=function(...) as.integer(0)

  if(is.null(nullInterval) | nullModel){
    doInterval = FALSE
    nullInterval = c(-1, 1)
    intervalCompl = FALSE
  }else{
    doInterval = TRUE
    intervalCompl = ifelse(!is.null(attr(nullInterval,"complement")),TRUE,FALSE)
    nullInterval = range(nullInterval)
  }

  r = cor(y, x, use="pairwise.complete.obs")
  n = length(x) - sum(is.na(y) | is.na(x))
  if(n != length(x)) message(paste("Ignored",sum(is.na(y) | is.na(x)),
                                "rows containing missing observations."))

  if(noSample){
    chains = matrix(as.numeric(NA),1,1)
  }else{
    if(nullModel) rscale = 0
    chains = metropCorrRcpp_jeffreys(r, n, 1/rscale, 1/rscale, TRUE,
                                     iterations, doInterval,
                                     .5*log((1+nullInterval)/(1-nullInterval)),
                                     intervalCompl, nullModel,
                                     progress, callback, callbackInterval)

    if(!nullModel & !noSample){
      acc.rate = mean(diff(chains) != 0)
      message("Independent-candidate M-H acceptance rate: ",round(100*acc.rate),"%")
    }
  }
  chains = mcmc(data.frame(chains))
  return(chains)
}


