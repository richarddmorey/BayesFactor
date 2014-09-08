
# constructor
BFBayesFactor <- function(numerator, denominator, bayesFactor, data){
  names(numerator) = rownames(bayesFactor)
  new("BFBayesFactor", numerator = numerator, 
      denominator = denominator, 
      bayesFactor = bayesFactor, 
      data = data, 
      version = BFInfo(FALSE))
}

setValidity("BFBayesFactor", function(object){
  if( length(object@numerator) != nrow(object@bayesFactor)) return("Number of numerator models does not equal number of Bayes factors.")
  numeratorsAreBFs = sapply(object@numerator,function(el) inherits(el,"BFmodel"))
  if( any(!numeratorsAreBFs)) return("Some numerators are not BFmodel objects.")
  # check numerators all have same data types as denominator
  dataTypeDenom = object@denominator@dataTypes
  dataTypesEqual = unlist(lapply(object@numerator, 
                function(model, compType) 
                  model@dataTypes %com% compType, 
                compType=dataTypeDenom))  
  if( any(!dataTypesEqual)) return("Data types are not equal across models.")
  
  typeDenom = object@denominator@type
  typesEqual = unlist(lapply(object@numerator, 
                             function(model, compType) 
                               identical(model@type, compType), 
                             compType=typeDenom))
  
  if( any(!typesEqual)) return("Model types are not equal across models.")

  classDenom = class(object@denominator)
  typesEqual = unlist(lapply(object@numerator, 
                             function(model, compType) 
                               identical(class(model), compType), 
                             compType=classDenom))
  
  if( any(!typesEqual)) return("Model classes are not equal across models.")
  
  # Check to see that Bayes factor data frame has required columns
  if( !all(colnames(object@bayesFactor) %in% c("bf", "error", "time", "code")) )
    return("Object does not have required columns (bf, error, time, code).")
  return(TRUE)
})


#' @rdname recompute-methods
#' @aliases recompute,BFBayesFactor-method
setMethod("recompute", "BFBayesFactor", function(x, progress = options()$BFprogress, ...){
  if(progress) lapply = pblapply 
  modelList = c(x@numerator,x@denominator)
  bfs = lapply(modelList, function(num, data, ...)
    compare(numerator = num, data = data, ...),
          data = x@data, progress = FALSE, ...)
  joined = do.call("c", bfs)
  numerators = joined[ 1:(length(joined)-1) ]
  denominator = joined[ length(joined) ]
  return(numerators / denominator)
})


#' @rdname BFBayesFactor-class
#' @name /,numeric,BFBayesFactor-method
#' @param e1 Numerator of the ratio
#' @param e2 Denominator of the ratio
setMethod('/', signature("numeric", "BFBayesFactor"), function(e1, e2){
  if( (e1 == 1) & (length(e2)==1) ){
    numer = e2@numerator[[1]]
    denom = list(e2@denominator)
    bf_df = e2@bayesFactor
    rownames(bf_df) = denom[[1]]@shortName
    bf_df$bf = -bf_df$bf
    bfobj = BFBayesFactor(numerator=denom, denominator=numer,
                          bayesFactor=bf_df, data=e2@data)
    return(bfobj)
  }else if( e1 != 1 ){
    stop("Dividend must be 1 (to take reciprocal).")
  }else if( length(e2)>1 ){
    allNum = as(e2,"list")
    BFlist = BFBayesFactorList(lapply(allNum, function(num) 1 / num))    
  }
}
)

#' @rdname BFBayesFactor-class
#' @name /,BFBayesFactor,BFBayesFactor-method
setMethod('/', signature("BFBayesFactor", "BFBayesFactor"), function(e1, e2){
    if( !(e1@denominator %same% e2@denominator) ) stop("Bayes factors have different denominator models; they cannot be compared.")
    if( !identical(e1@data, e2@data) ) stop("Bayes factors were computed using different data; they cannot be compared.")
    if( (length(e2)==1) ){
      errorEst = sqrt(e1@bayesFactor$error^2 + e2@bayesFactor$error^2)
      bfs = data.frame(bf=e1@bayesFactor$bf - e2@bayesFactor$bf,
                       error = errorEst, 
                       time  = date(), 
                       code = randomString(length(e1)))
      rownames(bfs) = rownames(e1@bayesFactor)
      
      # when bayes factors were computed at the same time or for the same model, 
      # they must be exactly equal  (no error)
      sameModel <- sapply(e1@numerator, function(num, den) num %same% den, den = e2@numerator[[1]])
      sameCode <- as.character(e1@bayesFactor$code) == as.character(e2@bayesFactor$code)
      bfs[ sameModel | sameCode, "error" ] = 0
      bfs[ sameModel | sameCode, "bf" ] = 0
      
      newbf = BFBayesFactor(numerator=e1@numerator,
                  denominator=e2@numerator[[1]],
                  bayesFactor=bfs,
                  data = e1@data)
      return(newbf)
    }else{
      allDenom = as(e2,"list")
      BFlist = BFBayesFactorList(lapply(allDenom, function(denom, num) num / denom, num = e1))
      return(BFlist)
    }
  }
)


setMethod('show', "BFBayesFactor", function(object){
  cat("Bayes factor analysis\n--------------\n")
  bfs = extractBF(object, logbf=TRUE)
  bfs$bf = sapply(bfs$bf, expString)
  indices = paste("[",1:nrow(bfs),"]",sep="")
  
  # pad model names
  nms = paste(indices,rownames(bfs),sep=" ")
  maxwidth = max(nchar(nms))
  nms = str_pad(nms,maxwidth,side="right",pad=" ")

  # pad Bayes factors
  maxwidth = max(nchar(bfs$bf))
  bfString = str_pad(bfs$bf,maxwidth,side="right",pad=" ")
  
  
  for(i in 1:nrow(bfs)){
    cat(nms[i]," : ",bfString[i]," \u00B1",round(bfs$error[i]*100,2),"%\n",sep="")
  }
  cat("\nAgainst denominator:\n")
  cat(" ",object@denominator@longName,"\n")
  cat("---\nBayes factor type: ",class(object@denominator)[1],", ",object@denominator@type,"\n\n",sep="")
  
})

setMethod('summary', "BFBayesFactor", function(object){
    show(object)
  })

#' @rdname BFBayesFactor-class
#' @name [,BFBayesFactor,index,missing,missing-method
#' @param x BFBayesFactor object
#' @param i indices indicating elements to extract
#' @param j unused for BFBayesFactor objects
#' @param drop unused
#' @param ... further arguments passed to related methods
setMethod("[", signature(x = "BFBayesFactor", i = "index", j = "missing",
                         drop = "missing"),
          function (x, i, j, ..., drop) {
            if((na <- nargs()) == 2){
              newbf = x
              x@numerator = x@numerator[i, drop=FALSE]            
              x@bayesFactor = x@bayesFactor[i, ,drop=FALSE]
            }else stop("invalid nargs()= ",na)
            return(x)
          })

#' @rdname extractBF-methods
#' @aliases extractBF,BFBayesFactor-method
setMethod("extractBF", "BFBayesFactor", function(x, logbf = FALSE, onlybf = FALSE){
  x = x@bayesFactor
  if(onlybf) x = x$bf
  if(!logbf) x$bf = exp(x$bf)
  return(x)
})


#' @rdname BFBayesFactor-class
#' @name t,BFBayesFactor-method
setMethod('t', "BFBayesFactor", function(x){
  return(t.BFBayesFactor(x))
})

setAs("BFBayesFactor", "data.frame",
      function( from, to ){
        as.data.frame.BFBayesFactor(from)
      })

setAs("BFBayesFactor" , "list",
       function ( from , to ){
          vec = vector(mode = "list", length = length(from) )
          for(i in 1:length(from))
            vec[i] = from[i]
          return(vec)
         })

setAs("BFBayesFactor", "vector",
      function( from, to ){
        as.vector.BFBayesFactor(from)
      })


#' @rdname BFBayesFactor-class
#' @name which.max,BFBayesFactor-method
setMethod("which.max", "BFBayesFactor", function(x)
  which.max.BFBayesFactor(x) )

#' @rdname BFBayesFactor-class
#' @name which.min,BFBayesFactor-method
setMethod("which.min", "BFBayesFactor", function(x)
  which.min.BFBayesFactor(x) )

#' @rdname BFBayesFactor-class
#' @name is.na,BFBayesFactor-method
setMethod("is.na", "BFBayesFactor", function(x)
  is.na.BFBayesFactor(x) )


######
# S3
######

##' This function coerces objects to the BFBayesFactor class
##' 
##' Function to coerce objects to the BFBayesFactor class
##' 
##' Currently, this function will only work with objects of class
##' \code{BFBayesFactorTop}, which are output from the functions \code{anovaBF}
##' and \code{regressionBF} when the \code{whichModels} argument is set to
##' \code{'top'}
##' @title Function to coerce objects to the BFBayesFactor class
##' @param object an object of appropriate class (for now, BFBayesFactorTop)
##' @return An object of class \code{BFBayesFactor}
##' @author Richard D. Morey (\email{richarddmorey@@gmail.com})
##' @export
##' @keywords misc
##' @seealso \code{\link{regressionBF}}, \code{anovaBF} whose output is 
##'   appropriate for use with this function when \code{whichModels='top'}
as.BFBayesFactor <- function(object)
  UseMethod("as.BFBayesFactor")


is.na.BFBayesFactor <- function(x){
  return(is.na(x@bayesFactor$bf))
}


names.BFBayesFactor <- function(x) {
  num <- sapply(x@numerator, function(el) el@shortName)
  den <- x@denominator@shortName
  return(list(numerator=num,denominator=den))  
}

length.BFBayesFactor <- function(x) 
  nrow(x@bayesFactor)

# See http://www-stat.stanford.edu/~jmc4/classInheritance.pdf
sort.BFBayesFactor <- function(x, decreasing = FALSE, ...){
  ord = order(x@bayesFactor$bf, decreasing = decreasing)
  return(x[ord])
}

max.BFBayesFactor <- function(..., na.rm=FALSE){
  joinedbf = do.call('c',list(...))
  el <- head(joinedbf, n=1)
  return(el)
}

min.BFBayesFactor <- function(..., na.rm=FALSE){
  joinedbf = do.call('c',list(...))
  el <- tail(joinedbf, n=1)
  return(el)
}

which.max.BFBayesFactor <- function(x){
  index = which.max(x@bayesFactor$bf)
  names(index) = rownames(x@bayesFactor)[index]
  return(index)
}

which.min.BFBayesFactor <- function(x){
  index = which.min(x@bayesFactor$bf)
  names(index) = rownames(x@bayesFactor)[index]
  return(index)
}

t.BFBayesFactor <- function(x){
  1/x
}


head.BFBayesFactor <- function(x, n=6L, ...){
  n = ifelse(n>length(x),length(x),n)
  x = sort(x, decreasing=TRUE)
  return(x[1:n])
}

tail.BFBayesFactor <- function(x, n=6L, ...){
  n = ifelse(n>length(x),length(x),n)
  x = sort(x)
  return(x[n:1])}

as.data.frame.BFBayesFactor <- function(x, row.names = NULL, optional=FALSE,...){
  df = x@bayesFactor
  df$bf = exp(df$bf) 
  return(df) 
}

as.vector.BFBayesFactor <- function(x, mode = "any"){
  if( !(mode %in% c("any", "numeric"))) stop("Cannot coerce to mode ", mode)  
  v = exp(x@bayesFactor$bf)
  names(v) = rownames(x@bayesFactor)
  return(v) 
}


c.BFBayesFactor <-
  function(..., recursive = FALSE)
  {
    z = list(...)
    if(length(z)==1) return(z[[1]])
    correctClass = unlist(lapply(z, function(object) inherits(object,"BFBayesFactor")))
    if(any(!correctClass)) stop("Cannot concatenate Bayes factor with non-Bayes factor.")
    
    dataAndDenoms = lapply(z, function(object){list(den = object@denominator, dat = object@data)})
    sameInfo = unlist(lapply(dataAndDenoms[-1], 
                             function(el, cmp){ 
                               sameType = el$dat %com% cmp$dat
                               sameType = ifelse(length(sameType)>0,sameType,TRUE)
                               (el$den %same% cmp$den) & sameType
                               },
                             cmp=dataAndDenoms[[1]]))
    if(any(!sameInfo)) stop("Cannot concatenate Bayes factors with different denominator models or data.")
    
    numerators = unlist(lapply(z, function(object){object@numerator}),recursive=FALSE, use.names=FALSE)    
    bfs = lapply(z, function(object){object@bayesFactor})
    df_rownames = unlist(lapply(z, function(object){rownames(object@bayesFactor)}))
    df_rownames = make.unique(df_rownames, sep=" #")
    bfs = do.call("rbind",bfs)
    rownames(bfs) = df_rownames
      
    bf = BFBayesFactor(numerator=numerators,
             denominator=z[[1]]@denominator, bayesFactor=bfs, 
             data = z[[1]]@data)
  }


