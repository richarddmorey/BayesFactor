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

BFproportion <- function(type, identifier, prior, shortName, longName){
  new("BFproportion", type = type,
      identifier = identifier,
      prior = prior,
      shortName = shortName,
      longName = longName,
      version = BFInfo(FALSE))
}


BFcontingencyTable <- function(type, identifier, prior, shortName, longName){
  new("BFcontingencyTable", type = type,
      identifier = identifier,
      prior = prior,
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

BFmetat <- function(type, identifier, prior, shortName, longName){
  new("BFmetat", type = type,
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

setMethod('show', signature = c("BFmodel"),
          function(object){
            cat("---\n Model:\n")
            cat("Type: ",class(object)[1],", ",object@type,"\n",sep="")
            cat(object@longName,"\n")
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
                complement = ifelse(!is.null(attr(nullInterval,"complement")),TRUE,FALSE)
                bf = ttest.tstat(t=t, n1=N,nullInterval=nullInterval,rscale=numerator@prior$rscale,complement=complement)
                numBF = bf[['bf']]
                errorEst = bf[['properror']]
              }
              numList = list(numerator)
              nms = numerator@shortName
              
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

setMethod('compare', signature(numerator = "BFindepSample", denominator = "missing", data = "data.frame"), 
          function(numerator, data, ...){
            
            formula = formula(numerator@identifier$formula)
            checkFormula(formula, data, analysis = "indept")
 
            dv = stringFromFormula(formula[[2]])            
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
                complement = ifelse(!is.null(attr(nullInterval,"complement")),TRUE,FALSE)
                bf = ttest.tstat(t=t, n1=ns[1], n2=ns[2], nullInterval=nullInterval,rscale=numerator@prior$rscale, complement = complement)
                numBF = bf[['bf']]
                errorEst = bf[['properror']]
              }
              numList = list(numerator)
              nms = numerator@shortName
      
              nullFormula = paste(formula[[2]],"1",sep=" ~ ")
              
              modDenominator = BFindepSample(type = "JZS", 
                                             identifier = list(formula = nullFormula), 
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
          

setMethod('compare', signature(numerator = "BFmetat", denominator = "missing", data = "data.frame"), 
          function(numerator, data, ...){
            
            nullInterval=numerator@prior$nullInterval

            if( (numerator@type=="JZS") ){
              if( numerator@identifier$formula=="d = 0" ){
                numBF = 0
                errorEst = 0
              }else{
                complement = ifelse(!is.null(attr(nullInterval,"complement")),TRUE,FALSE)
                bf = meta.ttest.tstat(t=data$t, n1=data$n1, n2=data$n2, 
                                    nullInterval=nullInterval, rscale=numerator@prior$rscale,complement=complement)
                numBF = bf[['bf']]
                errorEst = bf[['properror']]
              }
              numList = list(numerator)
              nms = numerator@shortName
              
              modDenominator = BFmetat(type = "JZS", 
                                           identifier = list(formula = "d = 0"), 
                                           prior=list(),
                                           shortName = "Null, d=0",
                                           longName = "Null, d = 0")
              
              bf_df = data.frame(bf = numBF,
                               error = errorEst,
                               time = date(),
                               code = randomString(length(numBF)))
              
            rownames(bf_df) <- nms
              
            newBF = BFBayesFactor(numerator = numList,
                                  denominator = modDenominator,
                                  data = data,
                                  bayesFactor = bf_df)
            return(newBF)
          }else{
            stop("Unknown prior type: ", numerator@type)
          }
})

setMethod('compare', signature(numerator = "BFproportion", denominator = "missing", data = "data.frame"), 
          function(numerator, data, ...){
            
            nullInterval=numerator@prior$nullInterval
            
            if( (numerator@type=="logistic") ){
              if( numerator@identifier$formula=="p = p0" ){
                numBF = 0
                errorEst = 0
              }else{
                complement = ifelse(!is.null(attr(nullInterval,"complement")),TRUE,FALSE)
                bf = prop.test.bf(y=data$y, N=data$N, p=numerator@prior$p0, 
                                  rscale=numerator@prior$rscale, nullInterval, 
                                  complement = complement)  
                numBF = bf[['bf']]
                errorEst = bf[['properror']]
              }
              numList = list(numerator)
              nms = numerator@shortName
              
              modDenominator = BFproportion(type = "logistic", 
                                       identifier = list(formula = "p = p0",p0=numerator@prior$p0), 
                                       prior=list(p0=numerator@prior$p0),
                                       shortName = paste("Null, p=",round(numerator@prior$p0,3),sep=""),
                                       longName = paste("Null, p = ", numerator@prior$p0, sep=""))
              
              bf_df = data.frame(bf = numBF,
                                 error = errorEst,
                                 time = date(),
                                 code = randomString(length(numBF)))
              
              rownames(bf_df) <- nms
              
              newBF = BFBayesFactor(numerator = numList,
                                    denominator = modDenominator,
                                    data = data,
                                    bayesFactor = bf_df)
              return(newBF)
            }else{
              stop("Unknown prior type: ", numerator@type)
            }
          })



setMethod('compare', signature(numerator = "BFcontingencyTable", denominator = "missing", data = "data.frame"), 
          function(numerator, data, ...){
            
            type = numerator@type
            a = numerator@prior$a
            marg = numerator@prior$fixedMargin
            data2 = as.matrix(data)
            if( !is.null(marg) )
              if( ( marg == "cols" ) & ( type == "independent multinomial" ) ) 
                data2 = t(data2)
              
            if(any(data%%1 != 0)) stop("All elements of x must be integers.")
            if(any(dim(data)<2) | (length(dim(data)) != 2)) stop("x must be m by n.")
            
            lbf = switch(type,
                         "poisson" = contingencyPoisson(as.matrix(data2), a),
                         "joint multinomial" = contingencyJointMultinomial(as.matrix(data2), a),
                         "independent multinomial" = contingencyIndepMultinomial(as.matrix(data2), a),
                         "hypergeometric" =  contingencyHypergeometric(as.matrix(data2), a),
                         stop("Unknown value of sampleType (see help for contingencyBF).")
            )
            error = 0
            
            denominator = BFcontingencyTable(type = type, 
                                             identifier = list(formula = "independence"), 
                                             prior=numerator@prior,
                                             shortName = paste0("Indep. (a=",a,")"),
                                             longName = paste0("Null, independence, a = ", a))
            
            bf_df = data.frame(bf = lbf,
                               error = error,
                               time = date(),
                               code = randomString(1))
            
            rownames(bf_df) <- numerator@shortName
            
            newBF = BFBayesFactor(numerator = list(numerator),
                                  denominator = denominator,
                                  data = as.data.frame(data),
                                  bayesFactor = bf_df
            )
            return(newBF)

})

setMethod('compare', signature(numerator = "BFcontingencyTable", denominator = "BFcontingencyTable", data = "data.frame"), 
          function(numerator, denominator, data, ...){
            if(!identical(numerator@type, denominator@type)) stop("Models of different types cannot be currently be compared by compare().")
            if(!identical(class(numerator), class(denominator))) stop("Models of different classes cannot be currently be compared by compare().")
                        
            BFnum = compare(numerator = numerator, data = data)              
            BFden = compare(numerator = denominator, data = data)
            return(BFnum / BFden)
          })

