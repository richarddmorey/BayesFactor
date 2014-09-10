proportionBF <- function(y, N, p, rscale = "medium", nullInterval = NULL, posterior=FALSE, callback = function(...) as.integer(0), ...)
{
  if(!is.null(nullInterval)){
    if(any(nullInterval<0) | any(nullInterval>1)) stop("nullInterval endpoints must be in [0,1].")
    nullInterval = range(nullInterval)
  }
  rscale = rpriorValues("proptest",,rscale)
  
  if( length(p) > 1 ) stop("Only a single null allowed (length(p) > 1).")
  if( length(y) != length(N) ) stop("Length of y and N must be the same.")
  if( any(y>N) | any(y < 0) ) stop("Invalid data (y>N or y<0).")
  if( any( c(y,N)%%1 != 0 ) ) stop("y and N must be integers.")
  
  hypNames = makePropHypothesisNames(rscale, nullInterval, p)
  
  mod1 = BFproportion(type = "logistic", 
                 identifier = list(formula = "p =/= p0", nullInterval = nullInterval, p0 = p), 
                 prior=list(rscale=rscale, nullInterval = nullInterval, p0 = p),
                 shortName = hypNames$shortName,
                 longName = hypNames$longName
  )      
  
  data = data.frame(y = y, N = N)
  
  if(posterior)
    return(posterior(mod1, data = data, callback = callback, ...))
  
  bf1 = compare(numerator = mod1, data = data)
  
  if(!is.null(nullInterval)){
    mod2 = mod1
    attr(mod2@identifier$nullInterval, "complement") = TRUE
    attr(mod2@prior$nullInterval, "complement") = TRUE
    hypNames = makePropHypothesisNames(rscale, mod2@identifier$nullInterval,p)
    mod2@shortName = hypNames$shortName
    mod2@longName = hypNames$longName
    
    bf2 = compare(numerator = mod2, data = data)
    return(c(bf1, bf2))
  }else{
    return(c(bf1))
  }   
}
