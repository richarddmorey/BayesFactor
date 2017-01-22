# constructor
BFodds <- function(BFinit, logodds = NULL, bayesFactor = NULL){
  if(is.null(logodds)) logodds = data.frame(odds = BFinit@bayesFactor$bf * 0)
  rownames(logodds) = rownames(BFinit@bayesFactor)
  new("BFodds", numerator = BFinit@numerator,
      denominator = BFinit@denominator,
      logodds = logodds,
      bayesFactor = bayesFactor,
      version = BFInfo(FALSE))
}

setValidity("BFodds", function(object){
  if( length(object@numerator) != nrow(object@logodds)) return("Number of numerator models does not equal number of Bayes factors.")
  if( !is.null(object@bayesFactor)){
    if( length(object@numerator) != length(object@bayesFactor)) return("Number of numerator models does not equal number of Bayes factors.")
    for(i in 1:length(object@numerator)){
      if( !(object@numerator[[i]] %same% object@bayesFactor@numerator[[i]])){
        return("Models in numerator are not the same as the numerators in Bayes factor.")
      }
    }
    if( !(object@denominator %same% object@denominator))
      return("Model in denominator is not the same as the denominator in Bayes factor.")
  }
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

  return(TRUE)
})


#' @rdname extractOdds-methods
#' @aliases extractOdds,BFodds-method
setMethod("extractOdds", "BFodds", function(x, logodds = FALSE, onlyodds = FALSE){
  z = x@logodds
  if(is.null(x@bayesFactor)){
    z$error = 0
  }else{
    bfs = extractBF(x@bayesFactor, logbf=TRUE)
    z$odds = z$odds + bfs$bf
    z$error = x@bayesFactor@bayesFactor$error
  }
  if(!logodds) z$odds = exp(z$odds)
  if(onlyodds) z = z$odds
  return(z)
})

setMethod('show', "BFodds", function(object){
  is.prior = is.null(object@bayesFactor)
  if(is.prior){
    cat("Prior odds\n--------------\n")
  }else{
    cat("Posterior odds\n--------------\n")
  }
  odds = extractOdds(object, logodds = TRUE)
  odds$odds = sapply(odds$odds, expString)

  indices = paste("[",1:nrow(odds),"]",sep="")

  # pad model names
  nms = paste(indices,rownames(odds),sep=" ")
  maxwidth = max(nchar(nms))
  nms = str_pad(nms,maxwidth,side="right",pad=" ")

  # pad Bayes factors
  maxwidth = max(nchar(odds$odds))
  oddsString = str_pad(odds$odds,maxwidth,side="right",pad=" ")



  for(i in 1:nrow(odds)){
    if(is.prior){
      cat(nms[i]," : ",oddsString[i],"\n",sep="")
    }else{
      cat(nms[i]," : ",oddsString[i]," \u00B1",round(odds$error[i]*100,2),"%\n",sep="")
    }
  }

  cat("\nAgainst denominator:\n")
  cat(" ",object@denominator@longName,"\n")
  cat("---\nModel type: ",class(object@denominator)[1],", ",object@denominator@type,"\n\n",sep="")

})

setMethod('summary', "BFodds", function(object){
  show(object)
})

#' @rdname BFodds-class
#' @name /,numeric,BFodds-method
#' @param e1 Numerator of the ratio
#' @param e2 Denominator of the ratio
setMethod('/', signature("numeric", "BFodds"), function(e1, e2){
  if( (e1 == 1) & (length(e2)==1) ){
    numer = e2@numerator[[1]]
    denom = list(e2@denominator)
    odds_df = e2@logodds
    if(is.null(e2@bayesFactor)){
      bf = NULL
    }else{
      bf = 1/e2@bayesFactor
    }
    rownames(odds_df) = denom[[1]]@shortName
    odds_df$odds = -odds_df$odds
    oddsobj = new("BFodds",numerator=denom, denominator=numer,
                  bayesFactor=bf, logodds = odds_df,
                  version = BFInfo(FALSE))
    return(oddsobj)
  }else if( e1 != 1 ){
    stop("Dividend must be 1 (to take reciprocal).")
  }else if( length(e2)>1 ){
    allNum = as(e2,"list")
    #BFlist = BFBayesFactorList(lapply(allNum, function(num) 1 / num))
    stop("Length of odds object must be ==1 to take reciprocal.")
  }
}
)

#' @rdname BFodds-class
#' @name /,BFodds,BFodds-method
setMethod('/', signature("BFodds", "BFodds"), function(e1, e2){
  if( length(e2) > 1) stop("Length of divisor must be ==1 to divide.")
  if( !(e1@denominator %same% e2@denominator) )
    stop("Odds have different denominator models; they cannot be compared.")
  if(!is.null(e1@bayesFactor) & !is.null(e1@bayesFactor)){
    bf = e1@bayesFactor / e2@bayesFactor
  }else if(is.null(e1@bayesFactor) & is.null(e1@bayesFactor)){
    bf = NULL
  }else{
    stop("Both odds objects must be prior, or both must be posterior.")
  }
  if( (length(e2)==1) ){
    logodds = data.frame(odds=e1@logodds$odds - e2@logodds$odds)
    rownames(logodds) = rownames(e1@logodds)

    oddsobj = new("BFodds",numerator=e1@numerator, denominator=e2@numerator[[1]],
                  bayesFactor=bf, logodds = logodds,
                  version = BFInfo(FALSE))

    return(oddsobj)
  }else{
    stop("Length of divisor must be ==1 to divide.")
  }
}
)

#' @rdname BFodds-class
#' @name *,BFodds,BFBayesFactor-method
setMethod('*', signature("BFodds", "BFBayesFactor"), function(e1, e2){
  if(!is.null(e1@bayesFactor))
    stop("Cannot multiply posterior odds object with Bayes factor.")
  new("BFodds", numerator = e1@numerator,
      denominator = e1@denominator,
      logodds = e1@logodds,
      bayesFactor = e2,
      version = BFInfo(FALSE))
}
)



#' @rdname BFodds-class
#' @name [,BFodds,index,missing,missing-method
#' @param x BFodds object
#' @param i indices indicating elements to extract
#' @param j unused for BFodds objects
#' @param drop unused
#' @param ... further arguments passed to related methods
setMethod("[", signature(x = "BFodds", i = "index", j = "missing",
                         drop = "missing"),
          function (x, i, j, ..., drop) {
            if((na <- nargs()) == 2){
              newodds = x
              x@numerator = x@numerator[i, drop=FALSE]
              x@logodds = x@logodds[i, ,drop=FALSE]
              if(is.null(x@bayesFactor)){
                x@bayesFactor = NULL
              }else{
                x@bayesFactor = x@bayesFactor[i]
              }
            }else stop("invalid nargs()= ",na)
            return(x)
          })

#' @rdname recompute-methods
#' @aliases recompute,BFodds-method
setMethod("recompute", "BFodds", function(x, progress = getOption('BFprogress', interactive()), multicore = FALSE, callback = function(...) as.integer(0), ...){
  bf = as.BFBayesFactor(x)
  bf = recompute(bf, progress = progress,
            multicore = multicore,
            callback = callback,
            ...)
  x@bayesFactor = bf
  return(x)
  })

#' @rdname priorOdds-method
#' @name priorOdds<-,BFodds,numeric-method
#' @docType methods
#' @exportMethod
setReplaceMethod("priorOdds", signature(object = "BFodds", value = "numeric"), definition = function (object, value) {
  priorLogodds(object) <- log(value)
  object
})

#' @rdname priorLogodds-method
#' @name priorLogodds<-,BFodds,numeric-method
#' @docType methods
#' @exportMethod
setReplaceMethod("priorLogodds", signature(object = "BFodds", value = "numeric"), definition = function (object, value) {
  object@logodds$odds <- value
  object
})


setAs("BFodds", "BFBayesFactor",
      function( from, to ){
        as.BFBayesFactor.BFodds(from)
      })

setAs("BFodds", "BFprobability",
      function( from, to ){
        as.BFprobability.BFodds(from)
      })


######
# S3
######


as.BFBayesFactor.BFodds <- function(object){
  if(!is.null(object@bayesFactor)){
    return(object@bayesFactor)
  }else{
    stop("Cannot convert prior odds to Bayes factor; no data has been given.")
  }
}

as.BFprobability.BFodds <- function(object, normalize = NULL, lognormalize = NULL){
  if(is.null(lognormalize) & is.null(normalize)){
    lognormalize = 0
  }else if(is.null(lognormalize) & !is.null(normalize)){
    lognormalize = log(normalize)
  }else if(!is.null(normalize)){
    stop("Cannot specify foth normalize and lognormalize.")
  }
  return(BFprobability(object, lognormalize))
}


length.BFodds <- function(x)
  nrow(x@logodds)

c.BFodds <-
  function(..., recursive = FALSE)
  {
    z = list(...)
    if(length(z)==1) return(z[[1]])
    correctClass = unlist(lapply(z, function(object) inherits(object,"BFodds")))
    if(any(!correctClass)) stop("Cannot concatenate odds with non-odds object.")

    denoms = lapply(z, function(object){ object@denominator })
    samedenom = unlist(lapply(denoms[-1],
                             function(el, cmp){
                               el %same% cmp
                             },
                             cmp=denoms[[1]]))
    if(any(!samedenom)) stop("Cannot concatenate odds objects with different denominator models.")

    logodds = lapply(z, function(object){object@logodds})
    df_rownames = unlist(lapply(z, function(object){rownames(object@logodds)}))
    df_rownames = make.unique(df_rownames, sep=" #")
    logodds = do.call("rbind",logodds)
    rownames(logodds) = df_rownames

    ### Grab the Bayes factors
    is.prior = 1:length(z) * NA
    for(i in 1:length(is.prior)){
      is.prior[i] = is.null(z[[i]]@bayesFactor)
    }
    if(all(!is.prior)){
      bfs = lapply(z, function(object){ object@bayesFactor })
      bfs = do.call("c", bfs)
      bf = BFodds(bfs,
                  logodds = logodds,
                  bayesFactor = bfs)
    }else if(all(is.prior)){
      numerators = unlist(lapply(z, function(object){object@numerator}),recursive=FALSE, use.names=FALSE)

      bf = new("BFodds", numerator = numerators, denominator = z[[1]]@denominator,
               logodds = logodds, bayesFactor = NULL, version = BFInfo(FALSE))
    }else{
      stop("Cannot concatenate prior odds with posterior odds.")
    }

    return(bf)
  }


names.BFodds <- function(x) {
  rownames(extractOdds(x))
}

# See http://www-stat.stanford.edu/~jmc4/classInheritance.pdf
sort.BFodds <- function(x, decreasing = FALSE, ...){
  ord = order(extractOdds(x, logodds=TRUE)$odds, decreasing = decreasing)
  return(x[ord])
}

max.BFodds <- function(..., na.rm=FALSE){
  joinedodds = do.call('c',list(...))
  el <- head(joinedodds, n=1)
  return(el)
}

min.BFodds <- function(..., na.rm=FALSE){
  joinedodds = do.call('c',list(...))
  el <- tail(joinedodds, n=1)
  return(el)
}

which.max.BFodds <- function(x){
  index = which.max(extractOdds(x, logodds=TRUE)$odds)
  names(index) = names(x)[index]
  return(index)
}

which.min.BFodds <- function(x){
  index = which.min(extractOdds(x, logodds=TRUE)$odds)
  names(index) = names(x)[index]
  return(index)
}

head.BFodds <- function(x, n=6L, ...){
  n = ifelse(n>length(x),length(x),n)
  x = sort(x, decreasing=TRUE)
  return(x[1:n])
}

tail.BFodds <- function(x, n=6L, ...){
  n = ifelse(n>length(x),length(x),n)
  x = sort(x)
  return(x[n:1])}

as.data.frame.BFodds <- function(x, row.names = NULL, optional=FALSE,...){
  df = extractOdds(x)
  return(df)
}

as.vector.BFodds <- function(x, mode = "any"){
  if( !(mode %in% c("any", "numeric"))) stop("Cannot coerce to mode ", mode)
  v = extractOdds(x)$odds
  names(v) = names(x)
  return(v)
}

