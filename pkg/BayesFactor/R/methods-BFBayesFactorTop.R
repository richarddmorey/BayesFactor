BFBayesFactorTop <- function(bf){
  if( class(bf@denominator) != "BFlinearModel" )
    stop("BFBayesFactorTopcan only be created from linear model objects.")

  len = sapply(bf@numerator, function(m){ 
    fmla = formula(m@identifier$formula)
    length(attr(terms(fmla),"term.labels"))
  })
  len_denom = length(attr(terms(formula(bf@denominator@identifier$formula)),"term.labels"))
  
  if( all( len < len_denom ) ){
    return(new("BFBayesFactorTop", bf))
  }else if( any( len > len_denom ) ){
    biggest = which(len == max(len))
    if(length(biggest) != 1) stop("Could not determine full model.")
    if(length(bf)==1)
      return(new("BFBayesFactorTop", 1/bf[biggest]))
    return(new("BFBayesFactorTop", bf[-biggest]/bf[biggest]))
  }else{
    stop("Could not determine full model.")
  }

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
            return(BFBayesFactorTop(x[i]))
          })


setMethod('summary', "BFBayesFactorTop", function(object){
  show(object)
})

#' @rdname recompute-methods
#' @aliases recompute,BFBayesFactor-method
setMethod("recompute", "BFBayesFactorTop", function(x, progress = options()$BFprogress, multicore = FALSE, callback = function(...) as.integer(0), ...){
  bf = recompute(as.BFBayesFactor(x), progress, multicore, callback, ...)
  BFBayesFactorTop(bf)
})  

setAs("BFBayesFactorTop", "BFBayesFactor",
      function( from, to ){
        as.BFBayesFactor.BFBayesFactorTop(from)
      })


######## S3

sort.BFBayesFactorTop <- function(x, decreasing = FALSE, ...){
  x = as.BFBayesFactor(x)
  x = sort(x,decreasing = decreasing, ...)
  return(BFBayesFactorTop(x))
}

as.BFBayesFactor.BFBayesFactorTop <- function(object){
  BFBayesFactor(numerator=object@numerator,
                denominator=object@denominator,
                bayesFactor = object@bayesFactor,
                data = object@data)
}

length.BFBayesFactorTop <- function(x) 
  length(as.BFBayesFactor(x))

