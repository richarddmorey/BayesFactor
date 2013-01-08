# constructors
BFmodel <- function(type, identifier, prior, dataTypes, shortName, longName){
  new("BFmodel", type = type,
      identifier = identifier,
      prior = prior,
      dataTypes = dataTypes,
      shortName = shortName,
      longName = longName,
      version = BFInfo(FALSE))
}

BFlinearModel <- function(type, identifier, prior, dataTypes, shortName, longName){
  new("BFlinearModel", type = type,
      identifier = identifier,
      prior = prior,
      dataTypes = dataTypes,
      shortName = shortName,
      longName = longName,
      version = BFInfo(FALSE))
}

BFoneSample <- function(type, identifier, prior, shortName, longName){
  new("BFoneSample", type = type,
      identifier = identifier,
      prior = prior,
      shortName = shortName,
      longName = longName,
      version = BFInfo(FALSE))
}

BFindepSample <- function(type, identifier, prior, shortName, longName){
  new("BFindepSample", type = type,
      identifier = identifier,
      prior = prior,
      shortName = shortName,
      longName = longName,
      version = BFInfo(FALSE))
}

#######

setMethod('show', signature = c("BFlinearModel"),
  function(object){
    cat("---\n Model:\n")
    cat("Type: ",class(object)[1],", ",object@type,"\n",sep="")
    cat(object@longName,"\n")
    cat("Data types:\n")
    lapply(names(object@dataTypes),function(el) cat(el,": ",object@dataTypes[el],"\n") )
    cat("\n\n")
  }
)

setMethod("%same%", signature = c(x="BFmodel",y="BFmodel"),
          function(x,y){
            classesSame = identical(class(x),class(y))
            dataTypeSame = x@dataTypes %com% y@dataTypes
            slotSame = sapply(slotNames(x), function(el,x,y) identical(slot(x,el),slot(y,el)),
                   x=x,y=y)
            slotSame["dataTypes"] = ifelse(length(dataTypeSame)>0,dataTypeSame, TRUE)
            # exclude version
            slotSame = slotSame[names(slotSame)!="version"]
            return(all(slotSame) & classesSame)
          })

#' @rdname compare-methods
#' @aliases compare,BFlinearModel,BFlinearModel,data.frame-method
setMethod('compare', signature(numerator = "BFlinearModel", denominator = "BFlinearModel", data = "data.frame"), 
          function(numerator, denominator, data, ...){
              if(!identical(numerator@type, denominator@type)) stop("Models of different types cannot be currently be compared by compare().")
              if(!identical(numerator@prior$mu, denominator@prior$mu)) stop("Models of different null means cannot currently be compared by compare()")
              if(!identical(class(numerator), class(denominator))) stop("Models of different classes cannot be currently be compared by compare().")
 
              LHSnum = all.vars(update(formula(numerator@identifier$formula), .~0))
              LHSden = all.vars(update(formula(denominator@identifier$formula), .~0))
              if(!identical(LHSnum, LHSden)) stop("Models have different dependent variables!")
              
              BFnum = compare(numerator = numerator, data = data)              
              BFden = compare(numerator = denominator, data = data)
              return(BFnum / BFden)
              })

#' @rdname compare-methods
#' @aliases compare,BFoneSample,missing,data.frame-method
setMethod('compare', signature(numerator = "BFoneSample", denominator = "missing", data = "data.frame"), 
          function(numerator, data, ...){

            formula = formula(numerator@identifier$formula)
            LHSnum = all.vars(update(formula, .~0))
            
            y = data[[LHSnum]]
            N = length(y)
            mu = numerator@prior$mu
            nullInterval=numerator@prior$nullInterval
            
            if( (numerator@type=="JZS") ){

              if( attr(terms(formula, data = data),"intercept") == 0 ){
                numBF = 0
                errorEst = 0
              }else{
                t = (mean(y) - mu) / sd(y) * sqrt(N)
                bf = ttest.tstat(t=t, n1=N,nullInterval=nullInterval,rscale=numerator@prior$rscale)
                numBF = bf[['bf']]
                errorEst = bf[['properror']]
              }
              if(!is.null(nullInterval)){

                modComplement = numerator
                modComplement@shortName = paste("Alt., r=",round(numerator@prior$rscale,3)," !(",nullInterval[1],"<d<",nullInterval[2],")",sep="")
                modComplement@longName = paste("Alternative, r = ",numerator@prior$rscale,", mu =/= ",mu, " !(",nullInterval[1],"<d<",nullInterval[2],")",sep="")
              
                numList = list(numerator,modComplement)
                nms = c(numerator@shortName,modComplement@shortName)
              }else{
                numList = list(numerator)
                nms = numerator@shortName
              }
              modDenominator = BFoneSample(type = "JZS", 
                                   identifier = list(formula = "y ~ 0"), 
                                   prior=list(mu=mu),
                                   shortName = paste("Null, mu=",mu,sep=""),
                                   longName = paste("Null, mu = ",mu, sep="")
              )
              
              bf_df = data.frame(bf = numBF,
                                 error = errorEst,
                                 time = date(),
                                 code = randomString(length(numBF)))
              
              rownames(bf_df) <- nms
              
              newBF = BFBayesFactor(numerator = numList,
                          denominator = modDenominator,
                          data = data,
                          bayesFactor = bf_df
              )
              return(newBF)
            }else{
              stop("Unknown prior type: ", numerator@type)
            }
          })

#' @rdname compare-methods
#' @aliases compare,BFindepSample,missing,data.frame-method
setMethod('compare', signature(numerator = "BFindepSample", denominator = "missing", data = "data.frame"), 
          function(numerator, data, ...){
            
            formula = formula(numerator@identifier$formula)
            checkFormula(formula, data, analysis = "indept")
 
            dv = deparse(formula[[2]])            
            factor = fmlaFactors(formula, data)[-1]
                        
            y = data[[dv]]
            iv = data[[factor]]
            ns = table(iv)
            
            mu = numerator@prior$mu
            nullInterval=numerator@prior$nullInterval
            
            if( mu != 0 ) stop("Indep. groups t test with nonzero null not supported yet.")
            
            if( (numerator@type=="JZS") ){
              if( length(attr(terms(formula, data = data),"term.labels")) == 0 ){
                numBF = 0
                errorEst = 0
              }else{
                t = t.test(formula = formula,data=data, var.eq=TRUE)$statistic
                bf = ttest.tstat(t=t, n1=ns[1], n2=ns[2], nullInterval=nullInterval,rscale=numerator@prior$rscale)
                numBF = bf[['bf']]
                errorEst = bf[['properror']]
              }
              if(!is.null(nullInterval)){
                
                modComplement = numerator
                modComplement@shortName = paste("Alt., r=",round(numerator@prior$rscale,3)," !(",nullInterval[1],"<d<",nullInterval[2],")",sep="")
                modComplement@longName = paste("Alternative, r = ",numerator@prior$rscale,", mu1-mu2 =/= ",mu, " !(",nullInterval[1],"<d<",nullInterval[2],")",sep="")
                
                numList = list(numerator,modComplement)
                nms = c(numerator@shortName,modComplement@shortName)
              }else{
                numList = list(numerator)
                nms = numerator@shortName
              }
              modDenominator = BFindepSample(type = "JZS", 
                                             identifier = list(formula = "y ~ 1"), 
                                             prior=list(mu=mu),
                                             shortName = paste("Null, mu1-mu2=",mu,sep=""),
                                             longName = paste("Null, mu1-mu2 = ",mu, sep="")
              )
              
              bf_df = data.frame(bf = numBF,
                                 error = errorEst,
                                 time = date(),
                                 code = randomString(length(numBF)))
              
              rownames(bf_df) <- nms
              
              newBF = BFBayesFactor(numerator = numList,
                                    denominator = modDenominator,
                                    data = data,
                                    bayesFactor = bf_df
              )
              return(newBF)
            }else{
              stop("Unknown prior type: ", numerator@type)
            }
            
          })
          



