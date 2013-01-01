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

#######

setGeneric("compare", function(numerator, denominator, data, ...) standardGeneric("compare"))
setMethod('compare', signature(numerator = "BFoneSample", denominator = "BFoneSample", data = "data.frame"), 
          function(numerator, denominator, data, ...){
              if(!identical(numerator@type, denominator@type)) stop("Models of different types cannot be currently be compared by compare().")
              if(!identical(numerator@prior$mu, denominator@prior$mu)) stop("Models of different null means cannot currently be compared by compare()")
              
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
            
            if( (numerator@type=="JZS one sample") ){

              if( attr(terms(formula),"intercept") == 0 ){
                numBF = 0
                errorEst = 0
              }else{
                t = (mean(y) - mu) / sd(y) * sqrt(N)
                bf = ttest.Quad(t=t, n1=N,nullInterval=nullInterval,rscale=numerator@prior$rscale, logbf=TRUE, error.est=TRUE)
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
              modDenominator = BFoneSample(type = "JZS one sample", 
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


############### Linear models
setMethod('compare', signature(numerator = "BFlinearModel", denominator = "BFlinearModel", data = "data.frame"), 
          function(numerator, denominator, data, ...){
            if(!identical(numerator@type, denominator@type)) stop("Models of different types cannot be currently be compared by compare().")
            if(!identical(numerator@prior$mu, denominator@prior$mu)) stop("Models of different null means cannot currently be compared by compare()")
            
            LHSnum = all.vars(update(formula(numerator@identifier$formula), .~0))
            LHSden = all.vars(update(formula(denominator@identifier$formula), .~0))
            if(!identical(LHSnum, LHSden)) stop("Models have different dependent variables!")
            
            BFnum = compare(numerator = numerator, data = data)              
            BFden = compare(numerator = denominator, data = data)
            return(BFnum / BFden)
          })


setMethod('compare', signature(numerator = "BFlinearModel", denominator = "missing", data = "data.frame"), 
          function(numerator, data, ...){
            
            formula = formula(numerator@identifier$formula)
            LHSnum = all.vars(update(formula, .~0))
            RHSnum = all.vars(update(formula, 0~.))
              
            
            if(length(RHSnum) > 1){
              #N-way ANOVA or multiple regression
              stop("Nway ANOVA and multiple linear regression not implemented.")
            }else if(length(RHSnum) == 1){
              iv = data[[RHSnum]]
              if(is.factor(iv)){
                ivNLvls = nlevels(iv)
                if(ivNLvls>2){
                  freqs <- table(iv)
                  if( all(freqs == freqs[1]) ){
                    # balanced one-way ANOVA    
                    stop("Balanced one-way ANOVA not implemented.")
                  }else{
                    # unbalanced one-way ANOVA
                    # (remember to check random/fixed)
                    # NwayAOV
                    stop("Unbalanced one-way ANOVA not implemented.")
                  }
                }else if(ivNLvls==2){
                  # indep. grp t test
                  
                  bf = compareIndept(numerator=numerator,data=data)                
                  return(bf)
                  
                }else{ # Nothing
                  stop('Not enough levels (<2) in independent variable.')
                } 
              }else{
                # simple linear regression
                stop("Simple linear regression not implemented.")
              }
            }
            
           })

compareIndept <- function(numerator,data) {
  formula = formula(numerator@identifier$formula)
  LHSnum = all.vars(update(formula, .~0))
  RHSnum = all.vars(update(formula, 0~.))
  
  y = data[[LHSnum]]
  iv = data[[RHSnum]]
  ns = table(iv)
  
  mu = numerator@prior$mu
  nullInterval=numerator@prior$nullInterval
  
  if( attr(terms(formula),"intercept") == 0 ) stop("Indep. groups t test without intercepts not supported yet.")
  if( mu != 0 ) stop("Indep. groups t test with nonzero null not supported yet.")
  
  if( (numerator@type=="JZS independent samples") ){
    if( length(attr(terms(formula),"term.labels")) == 0 ){
      numBF = 0
      errorEst = 0
    }else{
      t = t.test(formula = formula,data=data)$statistic
      bf = ttest.Quad(t=t, n1=ns[1], n2=ns[1], nullInterval=nullInterval,rscale=numerator@prior$rscale, logbf=TRUE, error.est=TRUE)
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
    modDenominator = BFlinearModel(type = "JZS independent samples", 
                                 identifier = list(formula = "y ~ 1"), 
                                 prior=list(mu=mu),
                                 dataTypes = numerator@dataTypes,
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
  
}


