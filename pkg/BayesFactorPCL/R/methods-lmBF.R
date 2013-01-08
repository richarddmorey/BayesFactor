lmBF <- function(formula, data, whichRandom = NULL, rscaleFixed="medium",
                 rscaleRandom=1, rscaleCont=1, posterior=FALSE, ...)
{    
  checkFormula(formula, data, analysis="lm")
  dataTypes <- createDataTypes(formula, whichRandom = whichRandom, data = data, analysis="lm")
  rscales = list(fixed=rscaleFixed, random=rscaleRandom, continuous=rscaleCont)
  
  numerator = BFlinearModel(type = "JZS", 
                        identifier = list(formula = deparse(formula)), 
                        prior=list(rscale=rscales),
                        dataTypes = dataTypes,
                        shortName = paste(deparse(formula[[3]]),sep=""),
                        longName = paste(deparse(formula),sep="")
      )
  
  if(posterior){
    chains = posterior(numerator, data = data, ...)
    return(chains)
  }else{
    bf = compare(numerator = numerator, data = data,...)
    return(bf)    
  }
}


############### Linear models

setMethod('compare', signature(numerator = "BFlinearModel", denominator = "missing", data = "data.frame"), 
  function(numerator, data, ...){

    rscaleFixed = rpriorValues("allNways","fixed",numerator@prior$rscale[['fixed']])
    rscaleRandom = rpriorValues("allNways","random",numerator@prior$rscale[['random']])
    rscaleCont = rpriorValues("regression",,numerator@prior$rscale[['continuous']])
                        
    formula = formula(numerator@identifier$formula)
    checkFormula(formula, data, analysis = "lm")
            
    factors = fmlaFactors(formula)[-1]
    nFactors = length(factors)
    dataTypes = numerator@dataTypes
    relevantDataTypes = dataTypes[names(dataTypes) %in% factors]
              
    dv = deparse(formula[[2]])
    if(numerator@type != "JZS") stop("Unknown model type.")
            
    denominator = BFlinearModel(type = "JZS", 
                    identifier = list(formula = paste(dv,"~ 1")), 
                    prior=list(),
                    dataTypes = dataTypes,
                    shortName = paste("Intercept only",sep=""),
                    longName = paste("Intercept only", sep="")
                  )
                     
    if( nFactors == 0 ){
      numerator = denominator
      bf = c(bf = 0, properror = 0)
    }else if(all(relevantDataTypes == "continuous")){
      ## Regression
      reg = summary(lm(formula,data=data))
      R2 = reg[[8]]
      N = nrow(data)
      p = length(factors)
      bf = linearReg.R2stat(N,p,R2,rscale=rscaleCont,logbf=TRUE, error.est=TRUE)
    }else if(all(relevantDataTypes != "continuous")){
      # ANOVA or t test
      freqs <- table(data[[factors[1]]])
      nLvls <- length(freqs)
      rscale = ifelse(dataTypes[factors[1]] == "fixed", rscaleFixed, rscaleRandom)              
      if( (nFactors > 1) | ( (nFactors == 1) & any(freqs!=freqs[1]))){ 
        # Nway ANOVA or unbalanced one-way ANOVA
        bf = nWayFormula(formula=formula, data = data, 
                dataTypes = dataTypes,
                rscaleFixed = rscaleFixed,
                rscaleRandom = rscaleRandom,
                gibbs = FALSE, ...)          
      }else if(nLvls>2){          
        # Balanced one-way
        Fstat = summary(aov(formula, data=data))[[1]]["F value"][1,] 
        J = length(freqs)
        N = freqs[1]
        bf = oneWayAOV.Fstat(Fstat, N, J, rscale, logbf=TRUE, error.est=TRUE)                
      }else if(nLvls==2){
        # independent groups t
        t = t.test(formula = formula,data=data, var.eq=TRUE)$statistic
        bf = ttest.tstat(t=t, n1=freqs[1], n2=freqs[2],rscale=rscale*sqrt(2), logbf=TRUE, error.est=TRUE)
      }else{ # Nothing
        stop("Too few levels in independent variable: ",factors[1])
      }
    }else{
      # GLM
      stop("GLM not implemented.")
    }
            
          
    nBF = length(bf[['bf']])
          
    bf_df = data.frame(bf = bf[['bf']],
                error = bf[['properror']],
                time = date(),
                code = randomString(nBF)
            )
          
    rownames(bf_df) <- numerator@shortName
          
    newBF = BFBayesFactor(numerator = list(numerator),
                denominator = denominator,
                data = data,
                bayesFactor = bf_df
            )
    return(newBF)
             
    }
)