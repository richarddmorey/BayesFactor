BFBayesFactorTop <- function(bf){
  nms = names(bf)[[1]]
  biggest = which.max(nchar(nms))
  new("BFBayesFactorTop", bf[-biggest]/bf[biggest])
}
  
setValidity("BFBayesFactorTop", function(object){
  if(class(object@denominator) != "BFlinearModel") 
    return("BFBayesFactorTop objects can only currently be created from BFlinearModel-type models.")
  
  omitted = lapply(object@numerator, whichOmitted, full = object@denominator)
  lens = sapply(omitted, length)
  
  if(any(lens != 1)) return("Not all numerators are formed by removing one term from the denominator.")
  
  return(TRUE)
})

setMethod('show', "BFBayesFactorTop", function(object){
  omitted = unlist(lapply(object@numerator, whichOmitted, full = object@denominator))
  cat("Bayes factor top-down analysis\n--------------\n")
  bfs = extractBF(object, logbf=TRUE)
  bfs$bf = sapply(bfs$bf, expString)
  indices = paste("[",1:nrow(bfs),"]",sep="")
  
  # pad model names
  nms = omitted
  maxwidth = max(nchar(nms))
  nms = str_pad(nms,maxwidth,side="right",pad=" ")
  
  # pad Bayes factors
  maxwidth = max(nchar(bfs$bf))
  bfString = str_pad(bfs$bf,maxwidth,side="right",pad=" ")
  
  cat("When effect is omitted from",object@denominator@shortName,", BF is...\n")
  for(i in 1:nrow(bfs)){
    cat(indices[i]," Omit ",nms[i]," : ",bfString[i]," \u00B1",round(bfs$error[i]*100,2),"%\n",sep="")
  }
  cat("\nAgainst denominator:\n")
  cat(" ",object@denominator@longName,"\n")
  cat("---\nBayes factor type: ",class(object@denominator)[1],", ",object@denominator@type,"\n\n",sep="")
})

#' @rdname BFBayesFactor-class
#' @name [,BFBayesFactorTop,index,missing,missing-method
setMethod("[", signature(x = "BFBayesFactorTop", i = "index", j = "missing",
                         drop = "missing"),
          function (x, i, j, ..., drop) {
            x = as(x,"BFBayesFactor")
            denom = (1/x[1])/(1/x[1]) 
            subset = "["(x, i=i,...)
            subset = c(subset, denom)
            return(BFBayesFactorTop(subset))
          })


setMethod('summary', "BFBayesFactorTop", function(object){
  show(object)
})

setAs("BFBayesFactorTop", "BFBayesFactor",
      function( from, to ){
        as.BFBayesFactor.BFBayesFactorTop(from)
      })


######## S3

sort.BFBayesFactorTop <- function(x, decreasing = FALSE, ...){
  x = as.BFBayesFactor(x)
  x = sort(x,decreasing = decreasing, ...)
  denom = (1/x[1])/(1/x[1]) 
  x = c(x, denom)
  return(BFBayesFactorTop(x))
}

as.BFBayesFactor.BFBayesFactorTop <- function(object){
  BFBayesFactor(numerator=object@numerator,
                denominator=object@denominator,
                bayesFactor = object@bayesFactor,
                data = object@data)
}


