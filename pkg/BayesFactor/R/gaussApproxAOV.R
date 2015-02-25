Qg <- function(q,sumSq,N,XtCnX,CnytCnX,rscale,gMap,gMapCounts,priorX=NULL,incCont=0,limit=TRUE)
{
  qLimits = options()$BFapproxLimits  
  zz = jzs_log_marginal_posterior_logg(q, sumSq, N, XtCnX, CnytCnX, rscale, gMap, gMapCounts, priorX, incCont, limit, qLimits, which = 0)
  return(zz[["d0g"]])
}

dQg <- function(q,sumSq,N,XtCnX,CnytCnX,rscale,gMap, gMapCounts, priorX=NULL, incCont=0)
{
  zz = jzs_log_marginal_posterior_logg(q, sumSq, N, XtCnX, CnytCnX, rscale, gMap, gMapCounts, priorX, incCont, FALSE, c(-Inf,Inf), which = 1)
  return(zz[['d1g']])
}


d2Qg <- function(q,sumSq,N,XtCnX,CnytCnX,rscale,gMap,gMapCounts,priorX=NULL,incCont=0)
{  
  zz = jzs_log_marginal_posterior_logg(q, sumSq, N, XtCnX, CnytCnX, rscale, gMap, gMapCounts, priorX, incCont, FALSE, c(-Inf, Inf), which = 2)
  return(zz[['d2g']])
}


hessianQg <- function(q,sumSq,N,XtCnX,CnytCnX,rscale,gMap,gMapCounts,priorX=NULL,incCont=0)
{
  diag(d2Qg(q,sumSq,N,XtCnX,CnytCnX,rscale,gMap,gMapCounts,priorX,incCont))
}

Qg_nlm <- function(q,sumSq,N,XtCnX,CnytCnX,rscale,gMap,gMapCounts,priorX=NULL,incCont=0)
{  
  zz = jzs_log_marginal_posterior_logg(q, sumSq, N, XtCnX, CnytCnX, rscale, gMap, gMapCounts, priorX, incCont, FALSE, c(-Inf,Inf), which = 2)
  
  res = -zz[['d0g']]
  attr(res, "gradient") <- -zz[['d1g']]
  attr(res, "hessian") <- -zz[['d2g']]
  return(res)
}

gaussianApproxAOV <- function(y,X,rscale,gMap,incCont=0)
{
  optMethod = options()$BFapproxOptimizer
  
  # dumb starting values
  qs = rscale * 0 
  
  if(!incCont){
    priorX = matrix(1,0,0)
  }else{
    priorX = matrix((t(X)%*%X)[1:incCont, 1:incCont],nrow=incCont)
  }
  
  N = length(y)
  Cny = matrix(y - mean(y), N)
  CnX = t(t(X) - colMeans(X))
  XtCnX = t(CnX) %*% CnX
  CnytCnX = t(Cny)%*%CnX
  sumSq = var(y) * (N-1) 
  gMapCounts = table(gMap)
  
  if(optMethod=="optim"){
    opt = optim(qs, Qg, gr = dQg,control=list(fnscale=-1),method="BFGS",sumSq=sumSq,N=N,XtCnX=XtCnX,CnytCnX=CnytCnX,rscale=rscale,gMap=gMap,gMapCounts=gMapCounts,priorX=priorX,incCont=incCont)
    if(opt$convergence) stop("Convergence not achieved in optim: ",opt$convergence)
    mu = opt$par
    val = opt$value
  }else if(optMethod=="nlm"){
    opt = nlm(Qg_nlm, qs, sumSq=sumSq,N=N,XtCnX=XtCnX,CnytCnX=CnytCnX,rscale=rscale,gMap=gMap,gMapCounts=gMapCounts, priorX=priorX,incCont=incCont, hessian=FALSE, check.analyticals=FALSE)
    if(opt$code>2) stop("Convergence not achieved in nlm: ",opt$code)
    val = -opt$minimum
    mu = opt$estimate
  }else{
    stop("unknown method in gaussianApproxAOV: ",optMethod)
  }
  
  hess = hessianQg(mu,sumSq=sumSq,N=N,XtCnX=XtCnX,CnytCnX=CnytCnX,rscale=rscale,gMap=gMap,gMapCounts=gMapCounts,priorX=priorX,incCont=incCont)  
  sig2 = -1/diag(hess)
  return(list(mu=mu,sig=sqrt(sig2),val=val))
}

laplaceAOV <- function(y,X,rscale,gMap,incCont=0)
{
  apx = gaussianApproxAOV(y,X,rscale,gMap,incCont)
  approxVal = sum(dnorm(apx$mu,apx$mu,apx$sig,log=TRUE))
  apx$val - approxVal
  
}
