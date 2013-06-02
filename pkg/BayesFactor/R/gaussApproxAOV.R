d2dinvgamma1 <- function(x,a,b)
  (a+1)/(x^2) - 2*b/(x^3)

ddinvgamma1 <- function(x,a,b)
  -(a+1)/x + b/(x^2)

dinvgamma1 <- function(x,a,b)
  a * log(b) - lgamma(a) - (a+1)*log(x) - b/x

Qg <- function(q,y,Xm,rscale,gMap,priorX=NULL,incCont=0,limit=TRUE)
{
  
  qLimits = options()$BFapproxLimits
  # Check to make sure values don't get too small or large, otherwise algorithm might fail
  if( (any(q< qLimits[1]) | any(q > qLimits[2])) & limit ) return(-Inf)
  
  
  g = exp(q)
  if(length(q) != length(rscale)) stop("length mismatch: q and r")
  N = nrow(Xm)
  p = ncol(Xm)
  
  yTilde = matrix(y - mean(y),N)
  XTilde = t(t(Xm) - colMeans(Xm))
  ybar = mean(y)
  
  if(ncol(Xm)==1){
    if(incCont) stop("Inappropriate use of Gaussian approximation with single column, continuous X.")
    Ginv= 1/g[gMap]
    Vg = sum(XTilde^2) + Ginv
    logDetVg = log(Vg)
    yXVXy = sum(y*XTilde)^2/Vg
  }else{
    Ginv= diag(1/g[gMap])
    if(incCont) Ginv[1:incCont,1:incCont] = priorX / g[gMap[1]]
    Vg = t(XTilde)%*%XTilde + Ginv
    logDetVg = determinant(Vg)$modulus
    yXVXy = t(yTilde)%*%XTilde%*%solve(Vg)%*%t(XTilde)%*%yTilde
  }
  
  ans = 
    -0.5 * sum(log(g[gMap])) +
    - 0.5 * logDetVg + 
    0.5*(N-1) * log(sum(y^2) - N * ybar^2) - 
    0.5*(N-1)*log(sum(yTilde^2) - yXVXy) + 
    sum(dinvgamma1(g, .5, rscale^2/2)) +
    sum(log(g))
  
  return(as.numeric(ans))
}

dQg <- function(q,y,Xm,rscale,gMap, priorX=NULL, incCont=0){
  g = exp(q)
  if(length(q) != length(rscale)) stop("length mismatch: q and r")
  N = dim(Xm)[1]
  nGs = max(gMap)
  yTilde = matrix(y - mean(y),N)
  XTilde = t(t(Xm) - colMeans(Xm))
  if(length(g)==1)
    stop("not implemented for single g")      
  
  Ginv= diag(1/g[gMap])
  if(incCont){
    Ginv[1:incCont,1:incCont] = priorX / g[gMap[1]]
  }
  
  ni = table(gMap)
  Vg = t(XTilde)%*%XTilde + Ginv
  Vginv = solve(Vg)

  yTildeXtilde = t(yTilde) %*% XTilde
  yTildeXtildeVginv = yTildeXtilde %*% Vginv
  num3 = sapply(1:nGs, function(index, gMap, yXVinv){
    sum(yXVinv[gMap==index]^2)
  }, gMap = gMap, yXVinv = yTildeXtildeVginv)
  
  tr2 = sapply(1:nGs, function(index, gMap, Vginv){
    sum(diag(Vginv)[gMap==index])
  }, gMap = gMap, Vginv = Vginv) 
  
  yXVXy = yTildeXtildeVginv %*% t(yTildeXtilde)
  
  ans = 
    -0.5 * ni / g +
    0.5 * tr2 / g^2 +       
    0.5*(N-1) * num3/(sum(yTilde^2) - yXVXy) / g^2 + 
    ddinvgamma1(g, .5, rscale^2/2) + 
    1/g
  
  return(ans*g)
  
}


d2Qg <- function(q,y,Xm,rscale,gMap,priorX=NULL,incCont=0){
  g = exp(q)
  if(length(q) != length(rscale)) stop("length mismatch: q and r")
  N = dim(Xm)[1]
  nGs = max(gMap)
  yTilde = matrix(y - mean(y),N)
  XTilde = t(t(Xm) - colMeans(Xm))
  if(length(g)==1)
    stop("not implemented for single g")
  ni = table(gMap)
  
  Ginv= diag(1/g[gMap])
  if(incCont){
    Ginv[1:incCont,1:incCont] = priorX / g[gMap[1]]
  }
  
  Vg = t(XTilde)%*%XTilde + Ginv
  Vginv = solve(Vg)
  Vginv2 = Vginv %*% Vginv
  
  yTildeXtilde = t(yTilde) %*% XTilde
  yTildeXtildeVginv = yTildeXtilde %*% Vginv
  yTildeXtildeVginv2 = yTildeXtilde %*% Vginv2
  
  num3 = sapply(1:nGs, function(index, gMap, yXVinv){
    sum(yXVinv[gMap==index]^2)
  }, gMap = gMap, yXVinv = yTildeXtildeVginv)
  
  num3.2 = sapply(1:nGs, function(index, gMap, yXVinv, yXVinv2){
    sum(yXVinv[gMap==index] * yXVinv2[gMap==index])
  }, gMap = gMap, yXVinv = yTildeXtildeVginv,yXVinv2 = yTildeXtildeVginv2)
  
  tr2 = sapply(1:nGs, function(index, gMap, Vginv){
    sum(diag(Vginv)[gMap==index])
  }, gMap = gMap, Vginv = Vginv)
  
  tr2.2 = sapply(1:nGs, function(index, gMap, Vginv2){
    sum(diag(Vginv2)[gMap==index])
  }, gMap = gMap, Vginv = Vginv2)
  
  yXVXy = yTildeXtildeVginv %*% t(yTildeXtilde)
  
  fg = num3
  dfg = 2 * num3.2 / g^2
  gg = 1/g^2
  dgg = -2/g^3
  hg = 1/(sum(yTilde^2) - yXVXy)
  dhg = num3 / g^2 * hg^2
  
  alpha = -0.5 * ni / g + 0.5 * tr2 / g^2 + 0.5*(N-1) * num3/(sum(yTilde^2) - yXVXy) / g^2 + ddinvgamma1(g, .5, rscale^2/2) + 1/g
  
  ans = 
    0.5 * ni / (g^2) +
    0.5 * (tr2.2 / g^4 - 2*tr2/g^3) +       
    0.5*(N-1) * (fg*gg*dhg + fg*hg*dgg + hg*gg*dfg) + 
    d2dinvgamma1(g, .5, rscale^2/2) +
    -1/g^2
  
  return(g * (ans*g + alpha))
  
}



hessianQg <- function(g,y,Xm,rscale,gMap,priorX=NULL,incCont=0){
  diag(d2Qg(g,y,Xm,rscale,gMap,priorX,incCont))
}

Qg_nlm <- function(q,y,Xm,rscale,gMap,priorX=NULL,incCont=0){
  g = exp(q)
  if(length(q) != length(rscale)) stop("length mismatch: q and r")
  N = dim(Xm)[1]
  nGs = max(gMap)
  yTilde = matrix(y - mean(y),N)
  XTilde = t(t(Xm) - colMeans(Xm))
  if(length(g)==1)
    stop("not implemented for single g")
  ni = table(gMap)
  
  Ginv= diag(1/g[gMap])
  if(incCont){
    Ginv[1:incCont,1:incCont] = priorX / g[gMap[1]]
  }
  
  
  Vg = t(XTilde)%*%XTilde + Ginv
  Vginv = solve(Vg)
  Vginv2 = Vginv %*% Vginv
  
  yTildeXtilde = t(yTilde) %*% XTilde
  yTildeXtildeVginv = yTildeXtilde %*% Vginv
  yTildeXtildeVginv2 = yTildeXtilde %*% Vginv2
  
  num3 = sapply(1:nGs, function(index, gMap, yXVinv){
    sum(yXVinv[gMap==index]^2)
  }, gMap = gMap, yXVinv = yTildeXtildeVginv)
  
  num3.2 = sapply(1:nGs, function(index, gMap, yXVinv, yXVinv2){
    sum(yXVinv[gMap==index] * yXVinv2[gMap==index])
  }, gMap = gMap, yXVinv = yTildeXtildeVginv,yXVinv2 = yTildeXtildeVginv2)
  
  tr2 = sapply(1:nGs, function(index, gMap, Vginv){
    sum(diag(Vginv)[gMap==index])
  }, gMap = gMap, Vginv = Vginv)
  
  tr2.2 = sapply(1:nGs, function(index, gMap, Vginv2){
    sum(diag(Vginv2)[gMap==index])
  }, gMap = gMap, Vginv = Vginv2)
  
  yXVXy = yTildeXtildeVginv %*% t(yTildeXtilde)
  
  ybar = mean(y)
  val = -0.5 * sum(log(g[gMap])) + -0.5 * determinant(Vg)$modulus + 
    0.5*(N-1) * log(sum(y^2) - N * ybar^2) - 
    0.5*(N-1)*log(sum(yTilde^2) - t(yTilde)%*%XTilde%*%solve(Vg)%*%t(XTilde)%*%yTilde) + 
    sum(dinvgamma1(g, .5, rscale^2/2)) +
    sum(q)
  
  dval = -0.5 * ni / g + 0.5 * tr2 / g^2 + 0.5*(N-1) * num3/(sum(yTilde^2) - yXVXy) / g^2 + ddinvgamma1(g, .5, rscale^2/2) + 1/g
  
  fg = num3
  dfg = 2 * num3.2 / g^2
  gg = 1/g^2
  dgg = -2/g^3
  hg = 1/(sum(yTilde^2) - yXVXy)
  dhg = num3 / g^2 * hg^2
  
  ddval = 
    0.5 * ni / (g^2) +
    0.5 * (tr2.2 / g^4 - 2*tr2/g^3) +       
    0.5*(N-1) * (fg*gg*dhg + fg*hg*dgg + hg*gg*dfg) + 
    d2dinvgamma1(g, .5, rscale^2/2) +
    -1/g^2
  
  grad = dval*g
  hess = diag(g^2 * ddval + grad)
  
  res = -as.numeric(val)
  attr(res, "gradient") <- -grad
  attr(res, "hessian") <- -hess
  return(res)
}


gaussianApproxAOV <- function(y,X,rscale,gMap,priorX=NULL,incCont=0){

  optMethod = options()$BFapproxOptimizer
  
  # gMap is written for C indexing
  gMap = gMap + 1
  # dumb starting values
  qs = rscale * 0 
  if(optMethod=="optim"){
    opt = optim(qs, Qg, gr = dQg,control=list(fnscale=-1),method="BFGS",y=y,Xm=X,rscale=rscale,gMap=gMap,priorX=priorX,incCont=incCont)
    if(opt$convergence) stop("Convergence not achieved in optim: ",opt$convergence)
    mu = opt$par
    val = opt$value
  }else if(optMethod=="nlm"){
    opt = nlm(Qg_nlm, qs, y=y,Xm=X,rscale=rscale,gMap=gMap, priorX=priorX,incCont=incCont, hessian=FALSE, check.analyticals=FALSE)
    if(opt$code>2) stop("Convergence not achieved in nlm: ",opt$code)
    val = -opt$minimum
    mu = opt$estimate
  }else{
    stop("unknown method in gaussianApproxAOV: ",optMethod)
  }
  
  hess = hessianQg(mu,y=y,Xm=X,rscale=rscale,gMap=gMap,priorX=priorX,incCont=incCont)  
  sig2 = -1/diag(hess)
  return(list(mu=mu,sig=sqrt(sig2),val=val))
}

laplaceAOV <- function(y,X,rscale,gMap,priorX=NULL,incCont=0)
{

  apx = gaussianApproxAOV(y,X,rscale,gMap,priorX,incCont)
  
  approxVal = sum(dnorm(apx$mu,apx$mu,apx$sig,log=TRUE))
  
  apx$val - approxVal
  
}

importanceAOV <- function(y,X,rscale,gMap,priorX=NULL,incCont=0,iterations = 10000)
{
  apx = gaussianApproxAOV(y,X,rscale,gMap,priorX,incCont)
  
  mu = apx$mu
  sig = apx$sig
  qSamp = replicate(iterations,{
    rnorm(mu,mu,sig)
  })

  sSamp = apply(qSamp,2,function(qs,y,Xm,rscale,gMap,mu,sig,priorX,incCont){
    exp(Qg(qs,y=y,Xm=Xm,rscale=rscale,gMap=gMap,priorX=priorX,incCont=incCont) - sum(dnorm(qs,mu,sig,log=TRUE)))
  },mu=mu,sig=sig,y=y,Xm=X,rscale=rscale,gMap=gMap,priorX=priorX,incCont=incCont)

  return(mean(sSamp))
}

