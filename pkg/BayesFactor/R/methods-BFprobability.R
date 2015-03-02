# constructor
BFprobability <- function(odds, normalize = 0){
  ## Add denominator 
  odds = c( odds, (1/odds[1]) / (1/odds[1]) )
  ## eliminate redundant models
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
  
  new("BFprobability", odds = odds, 
      normalize = normalize,
      version = BFInfo(FALSE))
}


setValidity("BFprobability", function(object){
  if( !is.numeric(object@normalize) )
    return("Normalization constant must be numeric.")
  odds = object@odds
  ## Add denominator 
  odds = c( odds, (1/odds[1]) / (1/odds[1]) )
  
  duplicates = 1:length(odds)
  for(i in 2:length(odds)){
    for(j in 1:(i-1)){
      if( odds@numerator[[i]] %same% odds@numerator[[j]] ){
        return("Duplicate models not allowed in probability objects.")
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
  odds = c(odds, (1/odds)/(1/odds))
  x = extractOdds(odds, log = TRUE)
  logsumodds = logMeanExpLogs(x$odds) + log(length(x$odds))
  logp = x$odds - logsumodds - norm
  z = data.frame(probs = logp, error = NA)
  rownames(z) = rownames(x)
  if(!logprobs) z$probs = exp(z$probs)
  if(onlyprobs) z = z$probs
  return(z)
})

######
# S3
######

length.BFprobability <- function(x) 
  length(x@odds) + 1

