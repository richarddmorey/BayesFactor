enumerateRegressionModels = function(fmla, whichModels, data){
  trms <- attr(terms(fmla, data = data), "term.labels")
  ntrms <- length(trms)
  dv = stringFromFormula(fmla[[2]])
  dv = composeTerm(dv)
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
  factors = composeTerms(factors)

  dv = stringFromFormula(formula[[2]])
  dv = composeTerm(dv)

  RHS = paste(factors,collapse=" + ")
  strng = paste(dv, " ~ ", RHS, collapse = "")
  return(formula(strng))
}

integrand.regression=Vectorize(function(g, N, p, R2, rscaleSqr=1, log=FALSE, log.const=0){
  a = .5 * ((N - p - 1 ) * log(1 + g) - (N - 1) * log(1 + g * (1 - R2)))
  shape=.5
  scale=rscaleSqr*N/2
  log.density.igam <- dinvgamma(g, shape, scale, log=TRUE)
  ans = a + log.density.igam - log.const
  ifelse(log,ans,exp(ans))
},"g")

# This is a more numerically stable version of the integrand, as a function of log(g)
integrand.regression.u=Vectorize(function(u, N, p, R2, rscaleSqr=1, log=FALSE, log.const=0, shift = 0){
  u = u + shift
  a = .5 * ((N - p - 1 ) * log1pExp(u) -
              (N - 1) * log1pExp(u + log(1 - R2)))
  shape=.5
  scale=rscaleSqr*N/2
  log.density.igam <- dinvgamma(u, shape, scale, log=TRUE, logx=TRUE)
  ans = a + log.density.igam - log.const + u
  ifelse(log,ans,exp(ans))
},"u")


linearReg.Gibbs <- function(y, covariates, iterations = 10000, rscale = "medium", progress = getOption('BFprogress', interactive()), callback=function(...) as.integer(0), noSample=FALSE, callbackInterval = 1, ignoreCols = NULL, thin = 1, ...){
  rscale = rpriorValues("regression",,rscale)
  X = apply(covariates,2,function(v) v - mean(v))
  y = matrix(y,ncol=1)
  N = length(y)

  P = ncol(X)
  nGs = 1

  if(is.null(ignoreCols)) ignoreCols = rep(0,P)

  # Check thinning to make sure number is reasonable
  if( (thin<1) | (thin>(iterations/3)) ) stop("MCMC thin parameter cannot be less than 1 or greater than iterations/3. Was:", thin)

  gMap = rep(0, P)

  sig2start = sum( (X%*%solve(t(X)%*%X)%*%t(X)%*%y - y)^2 ) / N

  progress = as.logical(progress)
  if(is.null(callback) | !is.function(callback)) callback=function(...) as.integer(0)

  if(noSample){ # Return structure of chains
    chains = matrix(NA,2,nOutputPars + 3 + nGs)
  }else{
    chains = jzs_Gibbs(iterations, y, cbind(1,X), rscale, 1, gMap, table(gMap), TRUE, FALSE,
                       as.integer(ignoreCols), as.integer(thin), as.logical(progress), callback, 1)
  }
  colnames(chains) = c("mu", colnames(covariates), "sig2", "g")
  chains = mcmc(chains)

  return(chains)

}

