if(getRversion() >= '2.15.1') globalVariables("gIndex")

mcoptions <- list(preschedule=FALSE, set.seed=TRUE)

combn2 <- function(x,lower=1){
  unlist(lapply(lower:length(x),function(m,x) combn(x,m,simplify=FALSE),x=x),recursive=FALSE)
}

fmlaFactors <- function(formula, data){
  rownames(attr(terms(formula, data = data),"factors"))
}

are.factors<-function(df) sapply(df, function(v) is.factor(v))

`%com%` <- function(x,y){
  common = intersect(names(x),names(y))
  if(length(common)==0) return(logical(0))
  all(sapply(common, function(el,x,y) identical(x[el],y[el]), x=x,y=y))
}

randomString <- function(x=1){
  n = ifelse(length(x)>1, length(x), x)
  substring(tempfile(rep("",n),"",""),2)
}

rpriorValues <- function(modelType,effectType=NULL,priorType=NULL){
  if(length(priorType)>1 | is.numeric(priorType)){
    return(priorType)
  }else if(length(priorType)==0){
    return(NULL)
  }
  
  if(modelType=="allNways"){
    return(
      switch(effectType,
             fixed = switch(priorType, 
                            wide=sqrt(2)/2, 
                            medium=1/2, 
                            stop("Unknown prior type.")),
             random = 1,
             stop("Unknown prior type.")
      )
    )
  }
  
  if(modelType=="ttest"){
    return(
      switch(priorType, wide=1, medium=sqrt(2)/2, stop("Unknown prior type."))  
    )
  }
  
  if(modelType=="regression"){
    return(1)  
  }
  
  stop("Unknown prior type.")
}


dinvgamma = function (x, shape, scale = 1) 
{
    if (shape <= 0 | scale <= 0) {
        stop("Shape or scale parameter negative in dinvgamma().\n")
    }
    alpha <- shape
    beta <- scale
    log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 
        1) * log(x) - (beta/x)
    return(exp(log.density))
}

# Taken from the WLE package source by Claudio Agostinelli <claudio at unive.it>
binary <- function(x, dim) {

   if (x==0) {
       pos <- 1
   } else {
       pos <- floor(log(x, 2))+1
   }

   if (!missing(dim)) {
       if (pos<=dim) {
           pos <- dim
       } else {
           warning("the value of `dim` is too small")
       }  
   }

   bin <- rep(0, pos)
   dicotomy <- rep(FALSE, pos)
   for (i in pos:1) {
        bin[i] <- floor(x/2^(i-1))
        dicotomy[i] <- bin[i]==1
        x <- x-((2^(i-1))*bin[i])
   }
   return(list(binary=bin, dicotomy=dicotomy))
}