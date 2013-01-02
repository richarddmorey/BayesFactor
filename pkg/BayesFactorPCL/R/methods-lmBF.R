lmBF <- function(formula, data, dataTypes, rscaleFixed="medium",
                 rscaleRandom=1,rscaleCont=1,...)
{    
  rscales = c(fixed=rscaleFixed, random=rscaleRandom, continuous=rscaleCont)
  
  numerator = BFlinearModel(type = "JZS", 
                        identifier = list(formula = deparse(formula)), 
                        prior=list(rscale=rscales),
                        dataTypes = dataTypes,
                        shortName = paste(deparse(formula[[3]]),sep=""),
                        longName = paste(deparse(formula),sep="")
      )
      
  bf = compare(numerator = numerator, data = data,...)
  return(bf)    
}


############### Linear models

setMethod('compare', signature(numerator = "BFlinearModel", denominator = "missing", data = "data.frame"), 
          function(numerator, data, ...){

            rscaleFixed = rpriorValues("allNways","fixed",numerator@prior$rscale['fixed'])
            rscaleRandom = rpriorValues("allNways","random",numerator@prior$rscale['random'])
            rscaleCont = rpriorValues("regression",,numerator@prior$rscale['continuous'])
            dataTypes = numerator@dataTypes
            dataRandom = getDataOfType("random",dataTypes,data) 
            dataFixed = getDataOfType("fixed",dataTypes,data)
            dataCont = getDataOfType("continuous",dataTypes,data)
            
            denominator = BFlinearModel(type = "JZS", 
                                         identifier = list(formula = "y ~ 1"), 
                                         prior=list(),
                                         dataTypes = dataTypes,
                                         shortName = paste("Intercept only",sep=""),
                                         longName = paste("Intercept only", sep="")
            )
            
            formula = formula(numerator@identifier$formula)
            if(attr(terms(formula),"response") == 0) stop("Dependent variable required in formula.")
            if(attr(terms(formula),"intercept") == 0) stop("Formula must include intercept.")
            
            nFactors = nrow(attr(terms(formula),"factors"))-1
            LHSnum = all.vars(update(formula, .~0))
            RHSnum = all.vars(update(formula, 0~.))
            allivs = data[RHSnum]
            if(numerator@type != "JZS") stop("Unknown model type.")
            
            if( nFactors == 0 ){
              numerator = denominator
              numBF = 0
              errorEst = 0
            }else if(nFactors > 1){
              #N-way ANOVA or multiple regression
              if(all(are.factors(allivs))){
                # N way aov
                bf = nWayAOV.MC(modelFormula=formula, y = data[[LHSnum]], 
                                dataFixed = dataFixed, 
                                dataRandom=dataRandom, 
                                rscaleFixed = rscaleFixed,
                                rscaleRandom = rscaleRandom,
                                logbf = TRUE, error.est= TRUE,
                                ...)
              }else  if(all(!are.factors(allivs))){
                reg = summary(lm(formula,data=data))
                R2 = reg[[8]]
                N = nrow(data)
                p = length(RHSnum)
                bf = linearReg.Quad(N,p,R2,rscale=rscaleCont,logbf=TRUE, error.est=TRUE)
              }else{
                stop("GLM not implemented.")
              }
            }else if(nFactors == 1){
              iv = data[[RHSnum]]
              if(is.factor(iv)){
                ivNLvls = nlevels(iv)
                freqs <- table(iv)
                if(ivNLvls>2){
                  if( all(freqs == freqs[1]) ){
                    # balanced one-way ANOVA    
                    Fstat = summary(aov(formula, data=data))[[1]]["F value"][1,] 
                    J = length(table(iv))
                    N = freqs[1]
                    if(numerator@dataTypes[[RHSnum]] == "fixed"){
                      bf = oneWayAOV.Quad(Fstat, N, J, rscaleFixed, logbf=TRUE, error.est=TRUE)
                    }else{
                      bf = oneWayAOV.Quad(Fstat, N, J, rscaleRandom, logbf=TRUE, error.est=TRUE)                      
                    }
                  }else{
                    # unbalanced one way AOV
                    bf = nWayAOV.MC(modelFormula=formula, y = data[[LHSnum]], 
                                    dataFixed = dataFixed, 
                                    dataRandom=dataRandom, 
                                    rscaleFixed = rscaleFixed,
                                    rscaleRandom = rscaleRandom,
                                    logbf = TRUE, error.est= TRUE,
                                    ...)
                  }
                }else if(ivNLvls==2){
                  t = t.test(formula = formula,data=data, var.eq=TRUE)$statistic
                  if(numerator@dataTypes[[RHSnum]] == "fixed"){
                    bf = ttest.Quad(t=t, n1=freqs[1], n2=freqs[2],rscale=rscaleFixed*sqrt(2), logbf=TRUE, error.est=TRUE) 
                  }else if(numerator@dataTypes[[RHSnum]] == "random"){
                    bf = ttest.Quad(t=t, n1=freqs[1], n2=freqs[2],rscale=rscaleRandom*sqrt(2), logbf=TRUE, error.est=TRUE) 
                  }else{
                    stop("Unknown data type.")  
                  }
                  
                }else{ # Nothing
                  stop('Not enough levels (<2) in independent variable.')
                } 
              }else{
                # simple linear regression
                reg = summary(lm(formula,data=data))
                R2 = reg[[8]]
                N = nrow(data)
                p = 1
                bf = linearReg.Quad(N,p,R2,rscale=rscaleCont,logbf=TRUE, error.est=TRUE)
                
              }
            }
            
          numBF = bf[['bf']]
          errorEst = bf[['properror']] 
          
          numList = list(numerator)
          nms = numerator@shortName
          
          bf_df = data.frame(bf = numBF,
                             error = errorEst,
                             time = date(),
                             code = randomString(length(numBF)))
          
          rownames(bf_df) <- nms
          
          newBF = BFBayesFactor(numerator = numList,
                                denominator = denominator,
                                data = data,
                                bayesFactor = bf_df
          )
          return(newBF)
          
            
            
          })