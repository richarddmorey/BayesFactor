if(getRversion() >= '2.15.1') globalVariables("gIndex")

mcoptions <- list(preschedule=FALSE, set.seed=TRUE)

# Create (new) factors out of factor and character columns
reFactorData <- function(data){
  if(is.data.frame(data)){
    indChar <- sapply(data, is.character)
    indFac <- sapply(data, is.factor)
    data[indChar | indFac] <- lapply(data[indChar | indFac], factor)
    return(data)
  }else{
    stop("Data must be in data.frame format.")
  }
}

filterVectorLogical <- function(columnFilter,myNames){
  if(!is.null(columnFilter)){
    ignoreMatrix = sapply(columnFilter, function(el,namedCols){
      grepl(el,namedCols)
    },namedCols=myNames)
    if(length(myNames)==1){
      ignoreCols = any(ignoreMatrix)
    }else{
      ignoreCols = apply(ignoreMatrix,1,any)
    }
    return(ignoreCols)
  }else{
    return(rep(FALSE,length(myNames)))
  }
}


expString <- function(x){
  if(is.na(x)) return("NA")
  doubleBase = .Machine$double.base
  toBase10log = x / log(10)
  toBaselog = x / log(doubleBase)

  numMax = .Machine$double.max.exp
  numMin = .Machine$double.min.exp

  if(toBaselog>numMax){
    first <- prettyNum( 10 ^ (toBase10log - floor(toBase10log)) )
    second <- prettyNum( floor(toBase10log) )
    return( paste( first, "e+", second, sep="" ) )
  }else if(toBaselog < numMin){
    first <- prettyNum( 10 ^ (1 - (ceiling(toBase10log) - toBase10log)) )
    second <- prettyNum( ceiling(toBase10log)-1 )
    return( paste( first, "e", second, sep="" ) )
  }else{
    return( prettyNum( exp(x) ) )
  }
}


alphabetizeTerms <- function(trms){
  splt = strsplit(trms,":",fixed=TRUE)
  sorted=lapply(splt, function(trm){
    if(length(trm)==1) return(trm)
    trm = sort(trm)
    paste(trm,collapse=":")
  })
  sorted = unlist(sorted)

  return(sorted)
}

whichOmitted <- function(numerator, full){
  fullFmla <- formula(full@identifier$formula)
  numFmla <- formula(numerator@identifier$formula)

  fullTrms <- attr(terms(fullFmla), "term.labels")
  numTrms <- attr(terms(numFmla), "term.labels")

  fullTrms = alphabetizeTerms(fullTrms)
  numTrms = alphabetizeTerms(numTrms)

  omitted = fullTrms[!(fullTrms %in% numTrms)]
  if(any( !(numTrms %in% fullTrms) )) stop("Numerator not a proper restriction of full.")
  return(omitted)
}


propErrorEst = function(logX){
  logX = logX[!is.na(logX)]
  summaries = logSummaryStats(logX)
  exp( ( summaries$logVar - log(length(logX)) )/2 - summaries$logMean)
}

combineModels <- function(modelList, checkCodes = TRUE){
  are.same = sapply(modelList[-1],function(m) modelList[[1]] %same% m)
  if( any(!are.same) ) stop("Cannot combine models that are not the same.")
  if(class(modelList[[1]]) != "BFlinearModel") return(modelList[[1]])

  hasanalysis = sapply(modelList, .hasSlot, name = "analysis")

  if( all(!hasanalysis) ) return(modelList[[1]])
  modelList = modelList[hasanalysis]
  if(length(modelList)==1) return(modelList[[1]])
  sampledTRUE = sapply(sapply(modelList, function(m) m@analysis[['sampled']]),identical,y=TRUE)


  if( !any(sampledTRUE)) return(modelList[[1]])
  modelList = modelList[which(sampledTRUE)]
  if( length(modelList)==1 ) return(modelList[[1]])

  bfs = unlist(sapply(modelList, function(m) m@analysis[['bf']]))
  properrs = unlist(sapply(modelList, function(m) m@analysis[['properror']]))

  # We need to make sure we don't combine analyses that are based on the same codes.
  codes = lapply(modelList, function(m) m@analysis[['code']])
  if(checkCodes){
    n = length(codes)
    X = diag(n)
    for(i in 2:n)
      for(j in 1:(i-1))
        X[i,j] = X[j,i] = length(intersect(codes[[i]],codes[[j]]))>0
    if(!identical(X,diag(n)))
      return(modelList[[which.min(properrs)]])
  }

  # Convert prop to abs err
  logAbs = bfs + log(properrs)
  # Compute log precisions
  logPrec = -2*logAbs
  # log sum of precisions
  logSumPrec = logMeanExpLogs(logPrec) + log(length(logPrec))
  # log weighted average
  logAvgBF = logMeanExpLogs(logPrec + bfs - logSumPrec) + log(length(logPrec))
  # convert prec back to abs err
  logSumAbs = -logSumPrec/2
  # convert back to prop err
  sumPropErr = exp(logSumAbs - logAvgBF)

  bf = logAvgBF
  properror = sumPropErr
  new.analysis = list(bf = bf, properror = properror, sampled = TRUE, method = "composite")

  all.codes = do.call("c",codes)

  new.mod = modelList[[1]]
  new.mod@analysis = new.analysis
  new.mod@analysis[['code']] = all.codes
  new.mod@version = BFInfo(FALSE)
  return(new.mod)
}

combn2 <- function(x,lower=1){
  unlist(lapply(lower:length(x),function(m,x) combn(x,m,simplify=FALSE),x=x),recursive=FALSE)
}

stringFromFormula <- function(formula){
  oneLine = paste(deparse(formula),collapse="")
  sub("\\s\\s+"," ", oneLine, perl=TRUE) # get rid of extra spaces
}

fmlaFactors <- function(formula, data){
  names <- rownames(attr(terms(formula, data = data),"factors"))
  names <- decomposeTerms(names)
  names <- unlist(names)
  names
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
  if(length(priorType)==0){
    return(NULL)
  }else if(length(priorType)>1 | is.numeric(priorType)){
    return(priorType)
  }else if( suppressWarnings( !is.na( as.numeric( priorType ) ) ) ){
    return(as.numeric(priorType))
  }else if(length(priorType)==0){
    return(NULL)
  }

  if(modelType=="proptest"){
    return(
      switch(priorType,
             ultrawide=1,
             wide=sqrt(2)/2,
             medium=1/2,
             stop("Unknown prior type."))
    )
  }

  if(modelType=="allNways"){
    return(
      switch(effectType,
             fixed = switch(priorType,
                            ultrawide=1,
                            wide=sqrt(2)/2,
                            medium=1/2,
                            stop("Unknown prior type.")),
             random = switch(priorType,
                             wide=sqrt(2)/2,
                             medium=1/2,
                             nuisance=1,
                             ultrawide=1,
                             stop("Unknown prior type.")),
             continuous = rpriorValues("regression",,priorType),
             stop("Unknown prior type.")
      )
    )
  }

  if(modelType=="ttestTwo"){
    return(
      switch(priorType,
             ultrawide=sqrt(2),
             wide=1,
             medium=sqrt(2)/2,
             stop("Unknown prior type."))
    )
  }

  if(modelType=="ttestOne"){
    return(
      switch(priorType,
             ultrawide=sqrt(2),
             wide=1,
             medium=sqrt(2)/2,
             stop("Unknown prior type."))
    )
  }


  if(modelType=="regression"){
    #return(1)
    return(
      switch(priorType,
             ultrawide=sqrt(2)/2,
             wide=1/2,
             medium=sqrt(2)/4,
             stop("Unknown prior type.")
      )
    )

  }

  if(modelType=="correlation"){
    return(
      switch(priorType,
        ultrawide=1,
        wide=1/sqrt(3),
        medium=1/3,
        medium.narrow = 1/sqrt(27),
        stop("Unknown prior type.")
        )
      )
   }


  stop("Unknown prior type.")
}


dinvgamma = function (x, shape, scale = 1, log = FALSE, logx = FALSE)
{
    if (shape <= 0 | scale <= 0) {
        stop("Shape or scale parameter negative in dinvgamma().\n")
    }
    shape = rep(0, length(x)) + shape
    scale = rep(0, length(x)) + scale
    if(logx){
      log.density = mapply(dinvgamma1_logx_Rcpp, x = x, a = shape, b = scale)
    }else{
      log.density = mapply(dinvgamma1_Rcpp, x = x, a = shape, b = scale)
    }
    if(log){
      return(log.density)
    }else{
      return(exp(log.density))
    }
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

# Construct all monotone Boolean functions for m arguments
monotoneBoolean <- function(m){
  if(m==0){
    return(list(FALSE,TRUE))
  }else{
    m0 = monotoneBoolean(m-1)
    m1 = list()
    for(i in 1:length(m0))
      for(j in 1:length(m0)){
        if(identical((m0[[i]] | m0[[j]]), m0[[j]])){
          m1[[length(m1)+1]] = c(m0[[i]],m0[[j]])
        }
      }
    return(m1)
  }
}

# Construct all monotone Boolean functions for m arguments
# but output in nice format (matrix)
monotoneBooleanNice = function(m){
  mb = monotoneBoolean(m)
  n = length(mb)
  mb = unlist(mb)
  dim(mb) = c(length(mb)/n,n)
  t(mb)
}

makeTerm <- function(m,factors){
  trms = factors[binary(m,length(factors))$dicotomy]
  trms = composeTerm(trms)
  trms
}

setMethod("%termin%", signature = c(x="character",table="character"),
          function(x,table){
            table = strsplit(table,":",fixed=TRUE)
            x = strsplit(x,":",fixed=TRUE)
            returnVector = rep(FALSE,length(x))
            for(i in 1:length(x))
              for(j in 1:length(table)){
                found = all(table[[j]] %in% x[[i]]) & all(x[[i]] %in% table[[j]])
                returnVector[i] = returnVector[i] | found
              }
            return(returnVector)
          })

setMethod("%termin%", signature = c(x="character",table="NULL"),
          function(x,table){
            return(rep(FALSE,length(x)))
           })


termMatch <- function(x, table, nomatch = NA_integer_){
  returnVector = rep(nomatch,length(x))
  if(is.null(table)){
    return(returnVector)
  }
  table = strsplit(table,":",fixed=TRUE)
  x = strsplit(x,":",fixed=TRUE)
  for(i in 1:length(x))
    for(j in 1:length(table)){
      found = all(table[[j]] %in% x[[i]]) & all(x[[i]] %in% table[[j]])
      if(is.na(returnVector[i]) & found) returnVector[i] = j
    }
  return(returnVector)
}

# Add two values for which the proportional error is known
# and return the proportional error
sumWithPropErr <- function(x1,x2,err1,err2){
  # convert proportional error to abs err
  logAbs1 = x1 + log(err1)
  logAbs2 = x2 + log(err2)
  logSum =  logExpXplusExpY( x1, x2 )
  absSum = .5 * logExpXplusExpY(2*logAbs1, 2*logAbs2)

  propErr = exp(absSum - logSum)
  return(c(logSum,propErr))
}

BFtry <- function(expression, silent=FALSE) {

  result <- base::try(expression, silent=silent)

  if (inherits(result, "try-error")) {

    message <- as.character(result)
    split <- base::strsplit(as.character(message), " : ")[[1]]
    error <- split[[length(split)]]

    while (substr(error, 1, 1) == ' ' || substr(error, 1, 1) == '\n')  # trim front
      error <- substring(error, 2)

    while (substring(error, nchar(error)) == ' ' || substring(error, nchar(error)) == '\n')  # trim back
      error <- substr(error, 1, nchar(error)-1)

    if (error == "Operation cancelled by callback function.")
      stop("Operation cancelled by callback function.")

    if (error == "Operation cancelled by interrupt.")
      stop("Operation cancelled by interrupt.")

  }

  result
}

marshallTibble <- function(data) {
    if (inherits(data, 'tbl_df')) {
        data <- as.data.frame(data)
        warning('data coerced from tibble to data frame', call.=FALSE)
    }
    data
}

# compose functions from jmvcore package

composeTerm <- function(components) {
  components <- sapply(components, function(component) {
    if (make.names(component) != component) {
      component <- gsub('\\', '\\\\', component, fixed=TRUE)
      component <- gsub('`', '\\`', component, fixed=TRUE)
      component <- paste0('`', component, '`')
    }
    component
  }, USE.NAMES=FALSE)
  term <- paste0(components, collapse=':')
  term
}

composeTerms <- function(listOfComponents) {
  sapply(listOfComponents, composeTerm, USE.NAMES=FALSE)
}

decomposeTerms <- function(terms) {
    decomposed <- list()
    for (i in seq_along(terms))
        decomposed[[i]] <- decomposeTerm(terms[[i]])
    decomposed
}

decomposeTerm <- function(term) {

    chars <- strsplit(term, '')[[1]]
    components <- character()
    componentChars <- character()
    inQuote <- FALSE

    i <- 1
    n <- length(chars)

    while (i <= n) {
        char <- chars[i]
        if (char == '`') {
            inQuote <- ! inQuote
        }
        else if (char == '\\') {
            i <- i + 1
            char <- chars[i]
            componentChars <- c(componentChars, char)
        }
        else if (char == ':' && inQuote == FALSE) {
            component <- paste0(componentChars, collapse='')
            components <- c(components, component)
            componentChars <- character()
        }
        else {
            componentChars <- c(componentChars, char)
        }
        i <- i + 1
    }

    component <- paste0(componentChars, collapse='')
    components <- c(components, component)

    components
}

