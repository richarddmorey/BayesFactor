# constructor
BFprobability <- function(odds, normalize = 0){
  ## Add denominator 
  
  if(options()$BFcheckProbabilityList){
    ## eliminate redundant models
    if( length(odds) > 1 ){
      odds = c( odds, (1/odds[1]) / (1/odds[1]) )
      duplicates = 1:length(odds)
      for(i in 2:length(odds)){
        for(j in 1:(i-1)){
          if( odds@numerator[[i]] %same% odds@numerator[[j]] ){
            duplicates[i] = j
            break
          }
        }
      }
      which.denom = duplicates[length(odds)]
      not.duplicate = duplicates == (1:length(odds))
      not.duplicate[ which.denom ] = FALSE
    
      # get rid of redundant models (this could be done better)  
      odds = odds[not.duplicate]
    }
  }
  new("BFprobability", odds = odds, 
      normalize = normalize,
      version = BFInfo(FALSE))
}


setValidity("BFprobability", function(object){
  if( !is.numeric(object@normalize) )
    return("Normalization constant must be numeric.")
  if( object@normalize > 0 )
    return("Normalization constant must be a valid probability.")
  odds = object@odds
  ## Add denominator 
  
  if(options()$BFcheckProbabilityList){
    if( length(odds) > 1 ){
      odds = c( odds, (1/odds[1]) / (1/odds[1]) )  
      duplicates = 1:length(odds)
      for(i in 2:length(odds)){
        for(j in 1:(i-1)){
          if( odds@numerator[[i]] %same% odds@numerator[[j]] ){
            return("Duplicate models not allowed in probability objects.")
          }
        }
      }
    }
  }
  return(TRUE)
})

setMethod('show', "BFprobability", function(object){
  odds = object@odds
  is.prior = is.null(object@odds@bayesFactor)
  if(is.prior){
    cat("Prior probabilities\n--------------\n")
  }else{
    cat("Posterior probabilities\n--------------\n")    
  }
  logprobs = extractProbabilities(object, logprobs = TRUE)
  logprobs$probs = sapply(logprobs$probs, expString)
  
  indices = paste("[",1:length(object),"]",sep="")
  
  # pad model names
  nms = paste(indices,rownames(logprobs),sep=" ")
  maxwidth = max(nchar(nms))
  nms = str_pad(nms,maxwidth,side="right",pad=" ")
  
  # pad Bayes factors
  maxwidth = max(nchar(logprobs$probs))
  probString = str_pad(logprobs$probs,maxwidth,side="right",pad=" ")
  
  for(i in 1:nrow(logprobs)){
    if(is.prior){
      cat(nms[i]," : ",probString[i],"\n",sep="")
    }else{
      cat(nms[i]," : ",probString[i]," \u00B1",round(logprobs$error[i]*100,2),"%\n",sep="")
    }
  }
  
  cat("\nNormalized probability: ", expString(object@normalize), " \n")
  cat("---\nModel type: ",class(object@odds@denominator)[1],", ",object@odds@denominator@type,"\n\n",sep="")
  
})

setMethod('summary', "BFprobability", function(object){
  show(object)
})

#' @rdname extractProbabilities-methods
#' @aliases extractProbabilities,BFprobability-method
setMethod("extractProbabilities", "BFprobability", function(x, logprobs = FALSE, onlyprobs = FALSE){
  norm = x@normalize
  odds = x@odds
  if( (length(odds) > 1 ) | !( odds@numerator[[1]] %same% odds@denominator ) ){
    odds = c(odds, (1/odds[1])/(1/odds[1]))
    x = extractOdds(odds, logodds = TRUE)
    logsumodds = logMeanExpLogs(x$odds) + log(length(x$odds))
    logp = x$odds - logsumodds + norm
    z = data.frame(probs = logp, error = NA)
  }else{ # numerator and denominator are the same
    x = extractOdds(odds, logodds = TRUE)
    z = data.frame(probs = norm, error = NA)
  }
  rownames(z) = rownames(x)
  if(!logprobs) z$probs = exp(z$probs)
  if(onlyprobs) z = z$probs
  return(z)
})

#' @rdname BFprobability-class
#' @name /,BFprobability,numeric-method
#' @param e1 BFprobability object
#' @param e2 new normalization constant
setMethod('/', signature("BFprobability", "numeric"), function(e1, e2){
  if(e2 > 1 | e2 <= 0)
    stop("Normalization constant must be >0 and not >1")
  return(e1 - log(e2))
}
)

#' @rdname BFprobability-class
#' @name -,BFprobability,numeric-method
setMethod('-', signature("BFprobability", "numeric"), function(e1, e2){
  if(length(e2)>1) stop("Normalization constant must be a scalar.")
  if(e2 > 0 | e2 == -Inf)
    stop("Normalization constant must be >0 and not >1")
  e1@normalize = e2
  return(e1)
}
)

#' @rdname BFprobability-class
#' @name [,BFprobability,index,missing,missing-method
#' @param x BFprobability object
#' @param i indices indicating elements to extract
#' @param j unused for BFprobability objects
#' @param drop unused
#' @param ... further arguments passed to related methods
setMethod("[", signature(x = "BFprobability", i = "index", j = "missing",
                         drop = "missing"),
          function (x, i, j, ..., drop) {
            if((na <- nargs()) == 2){
              if(is.logical(i)){
                if(any(i)){
                  i = (1:length(i))[i]
                }else{
                  return(NULL)
                }
              }
              i = unique(i)
              norm = x@normalize
              logprobs = extractProbabilities(x, logprobs = TRUE)[i, ,drop=FALSE]
              sumlogprobs = logMeanExpLogs(logprobs$probs) + log(nrow(logprobs))
              if(length(x) == length(i) ){
                newnorm = norm
              }else if( length(i) == 1){
                newnorm = sumlogprobs
              }else{
                newnorm = norm + sumlogprobs
              }
              whichnum = i[1:max(1, length(i)-1)]
              whichdenom = i[length(i)]
              newodds = c(x@odds, (1/x@odds[1])/(1/x@odds[1]))
              newodds = newodds[whichnum] / newodds[whichdenom]
              x = BFprobability( newodds,  newnorm )
            }else stop("invalid nargs()= ",na)
            return(x)
          })

#' @rdname BFprobability-class
#' @name filterBF,BFprobability,character-method
#' @param name regular expression to search name
#' @param perl logical. Should perl-compatible regexps be used? See ?grepl for details.
#' @param fixed logical. If TRUE, pattern is a string to be matched as is. See ?grepl for details.
setMethod("filterBF", signature(x = "BFprobability", name = "character"),
          function (x, name, perl, fixed, ...) {
            my.names = names(x) 
            matches = sapply(name, function(el){
              grepl(el, my.names, fixed = fixed, perl = perl)
            })
            any.matches = apply(matches, 1, any)
            x[any.matches]
          }
)
            


######
# S3
######

##' This function coerces objects to the BFprobability class
##' 
##' Function to coerce objects to the BFprobability class
##' 
##' Currently, this function will only work with objects of class
##' \code{BFOdds}.
##' @title Function to coerce objects to the BFprobability class
##' @param object an object of appropriate class (BFodds)
##' @param normalize the sum of the probabilities for all models in the object (1 by default)
##' @param lognormalize alternative to \code{normalize}; the 
##' logarithm of the normalization constant (0 by default)
##' @return An object of class \code{BFprobability}
##' @author Richard D. Morey (\email{richarddmorey@@gmail.com})
##' @export
##' @keywords misc
as.BFprobability <- function(object, normalize = NULL, lognormalize = NULL)
  UseMethod("as.BFprobability")


length.BFprobability <- function(x) 
  nrow(extractProbabilities(x))

names.BFprobability <- function(x) {
  rownames(extractProbabilities(x))
}

# See http://www-stat.stanford.edu/~jmc4/classInheritance.pdf
sort.BFprobability <- function(x, decreasing = FALSE, ...){
  ord = order(extractProbabilities(x, logprobs=TRUE)$probs, decreasing = decreasing)
  return(x[ord])
}

max.BFprobability <- function(..., na.rm=FALSE){
  if(nargs()>2) stop("Cannot concatenate probability objects.")
  el <- head(list(...)[[1]], n=1)
  return(el)
}

min.BFprobability <- function(..., na.rm=FALSE){
  if(nargs()>2) stop("Cannot concatenate probability objects.")
  el <- tail(list(...)[[1]], n=1)
  return(el)
}

which.max.BFprobability <- function(x){
  index = which.max(extractProbabilities(x, logprobs=TRUE)$probs)
  names(index) = names(x)[index]
  return(index)
}

which.min.BFprobability <- function(x){
  index = which.min(extractProbabilities(x, logprobs=TRUE)$probs)
  names(index) = names(x)[index]
  return(index)
  
}

head.BFprobability <- function(x, n=6L, ...){
  n = ifelse(n>length(x),length(x),n)
  x = sort(x, decreasing=TRUE)
  return(x[1:n])
}

tail.BFprobability <- function(x, n=6L, ...){
  n = ifelse(n>length(x),length(x),n)
  x = sort(x)
  return(x[n:1])}

as.data.frame.BFprobability <- function(x, row.names = NULL, optional=FALSE,...){
  df = extractProbabilities(x)
  return(df) 
}

as.vector.BFprobability <- function(x, mode = "any"){
  if( !(mode %in% c("any", "numeric"))) stop("Cannot coerce to mode ", mode)  
  v = extractProbabilities(x)$probs
  names(v) = names(x)
  return(v) 
}

sum.BFprobability <-
  function(..., na.rm = FALSE)
  {
    if(na.rm) warning("na.rm argument not used for BFprobability objects.")
    sapply(list(...), function(el){
      if(is(el, "BFprobability")){
        return(exp(el@normalize))
      }else{
        return(NA)
      }
    }, USE.NAMES = FALSE)
  }

