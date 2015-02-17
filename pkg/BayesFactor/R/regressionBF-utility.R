enumerateRegressionModels = function(fmla, whichModels, data){
  trms <- attr(terms(fmla, data = data), "term.labels")
  ntrms <- length(trms)
  dv = stringFromFormula(fmla[[2]])
  if(ntrms == 1 ) whichModels = "all"
  
  if(whichModels=="top"){
    lst = combn2( trms, ntrms - 1 )
  }else if(whichModels=="all"){
    lst = combn2( trms, 1 )
  }else if(whichModels=='bottom'){
    lst = as.list(combn( trms, 1 ))
  }else{
    stop("Unknown whichModels value: ",whichModels)
  }
  strng <- lapply(lst,function(el){
    paste(el,collapse=" + ")
  })
  fmla <- lapply(strng, function(el){
    formula(paste(dv,"~", el))
  })
  return(fmla)
}

createFullRegressionModel <- function(formula, data){
  factors = fmlaFactors(formula, data)[-1]
  
  dv = stringFromFormula(formula[[2]])
  
  RHS = paste(factors,collapse=" + ")
  strng = paste(dv, " ~ ", RHS, collapse = "")
  return(formula(strng))
}

integrand.regression=Vectorize(function(g, N, p , R2, rscaleSqr=1, log=FALSE, log.const=0){
  a = .5 * ((N - p - 1 ) * log(1 + g) - (N - 1) * log(1 + g * (1 - R2)))
  shape=.5
  scale=rscaleSqr*N/2
  log.density.igam <- dinvgamma1(g, shape, scale)
  ans = a + log.density.igam - log.const
  ifelse(log,ans,exp(ans))
},"g")

linearReg.Gibbs <- function(y, covariates, iterations = 10000, rscale = "medium", progress = options()$BFprogress, callback=function(...) as.integer(0), noSample=FALSE, callbackInterval = 1, ...){
  rscale = rpriorValues("regression",,rscale)
  X <- apply(covariates,2,function(v) v - mean(v))
  y = matrix(y,ncol=1)
  
  sig2start = sum( (X%*%solve(t(X)%*%X)%*%t(X)%*%y - y)^2 ) / N
  
  progress = as.logical(progress)
  if(is.null(callback) | !is.function(callback)) callback=function(...) as.integer(0)
  
  
  if(noSample){
    chains = matrix(NA,ncol(covariates)+2,2)
  }else{
    chains = GibbsLinearRegRcpp(as.integer(iterations), 
                       y, X, rscale, sig2start, FALSE, 
                       progress, callback, callbackInterval)
  }
  
  colnames(chains) = c(colnames(covariates),"sig2","g")
  return(mcmc(chains))
  
}

