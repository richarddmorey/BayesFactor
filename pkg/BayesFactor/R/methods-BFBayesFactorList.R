
BFBayesFactorList<- function(li){
  col_nms = sapply(li,function(el) el@denominator@shortName)
  names(li) = make.unique(col_nms, sep=" #")
  new("BFBayesFactorList", li, version=BFInfo(FALSE))
}

setValidity("BFBayesFactorList", function(object){
  firstNumerator = object[[1]]@numerator
  sameNumerators = unlist(lapply(object, function(el, firstNumerator) {
      identical(el@numerator,firstNumerator)
    }, firstNumerator = firstNumerator))
  if(any(!sameNumerators)) return("All numerators in elements of BayesFactorList must be identical")
  return(TRUE)
})

setMethod('show', "BFBayesFactorList", function(object){
  print(as(object,"matrix"))
})

#' @rdname BFBayesFactorList-class
#' @name t,BFBayesFactorList-method
#' @param x a BFBayesFactorList object
setMethod('t', "BFBayesFactorList", function(x){
  return(1/x)
})

#' @rdname BFBayesFactorList-class
#' @name /,numeric,BFBayesFactorList-method
#' @param e1 Numerator of the ratio
#' @param e2 Denominator of the ratio
setMethod('/', signature("numeric", "BFBayesFactorList"), function(e1, e2){
  if( (e1 == 1) & (length(e2[[1]])==1) ){
    bflist = lapply(e2,function(el) 1/el)
    return(do.call('c',bflist))
  }else if( e1 != 1 ){
    stop("Dividend must be 1 (to take reciprocal).")
  }else if( length(e2[[1]])>1 ){
    vec = vector(mode = "list", length = length(e2[[1]]))
    for(i in 1:length(e2[[1]])){
      vec[[i]] = 1/e2[i,]
    }
    bflist = BFBayesFactorList(vec)
    return(bflist)
  }
}
)

#' @rdname BFBayesFactorList-class
#' @name [,BFBayesFactorList,index,index,missing-method
#' @param i indices specifying rows to extract
#' @param j indices specifying columns to extract
#' @param drop unused
#' @param ... further arguments passed to related methods
setMethod("[", signature(x = "BFBayesFactorList", i = "index", j = "index",
                         drop = "missing"),
          function (x, i, j, ..., drop) {
            if((na <- nargs()) == 3){
              x = x[i,][,j]
            }else stop("invalid nargs()= ",na)
            return(x)
          })

#' @rdname BFBayesFactorList-class
#' @name [,BFBayesFactorList,index,missing,missing-method
setMethod("[", signature(x = "BFBayesFactorList", i = "index", j = "missing",
                         drop = "missing"),
          function (x, i, j, ..., drop) {
            if((na <- nargs()) == 3){
              bfs = lapply(x,function(el,i) el[i], i = i)
              x = BFBayesFactorList(bfs)
            }else stop("invalid nargs()= ",na)
            return(x)
          })

#' @rdname BFBayesFactorList-class
#' @name [,BFBayesFactorList,missing,index,missing-method
setMethod("[", signature(x = "BFBayesFactorList", i = "missing", j = "index",
                         drop = "missing"),
          function (x, i, j, ..., drop) {
            if((na <- nargs()) == 3){
              if(length(j)==1){
                x = x[[j]]
              }else if(length(j)>1){
                x = as(x, "vector")
                x = BFBayesFactorList(x[j])
              }
            }else stop("invalid nargs()= ",na)
            return(x)
          })


setAs("BFBayesFactorList" , "list",
      function ( from , to ){
        as.vector(from)
      })

setAs("BFBayesFactorList" , "vector",
      function ( from , to ){
        as.vector(from)
      })

setAs("BFBayesFactorList" , "matrix",
      function ( from , to ){
        as.matrix(from)
      })

## S3 Methods
#####

as.vector.BFBayesFactorList <- function(x, mode = "any"){
  if( !(mode %in% c("any", "list"))) stop("Cannot coerce to mode ", mode)
  vec = vector(mode = "list", length = length(x) )
  for(i in 1:length(x))
    vec[[i]] = x[[i]]
  names(vec) = names(x)
  return(vec)
}

as.matrix.BFBayesFactorList <- function(x,...){
  matr <- sapply(x, as.vector)
  dim(matr) <- c(length(x[[1]]),length(x))
  numNames <- rownames(extractBF(x[[1]]))
  denNames <- names(x)
  dimnames(matr) = list(numerator=numNames, denominator=denNames)
  return(as.matrix(matr))
}

