

############### Linear models

setMethod('compare', signature(numerator = "BFlinearModel", denominator = "missing", data = "data.frame"),
  function(numerator, data, ...){

    if(!.hasSlot(numerator,"analysis")) numerator@analysis = list()
    old.numerator = numerator
    rscaleFixed = rpriorValues("allNways","fixed",numerator@prior$rscale[['fixed']])
    rscaleRandom = rpriorValues("allNways","random",numerator@prior$rscale[['random']])
    rscaleCont = rpriorValues("regression",,numerator@prior$rscale[['continuous']])
    rscaleEffects = numerator@prior$rscale[['effects']]


    formula = formula(numerator@identifier$formula)
    checkFormula(formula, data, analysis = "lm")

    factors = fmlaFactors(formula, data)[-1]
    nFactors = length(factors)
    dataTypes = numerator@dataTypes
    relevantDataTypes = dataTypes[names(dataTypes) %in% factors]

    dv = stringFromFormula(formula[[2]])
    if(numerator@type != "JZS") stop("Unknown model type.")

    denominator = BFlinearModel(type = "JZS",
                    identifier = list(formula = paste(dv,"~ 1")),
                    prior=list(),
                    dataTypes = dataTypes,
                    shortName = paste("Intercept only",sep=""),
                    longName = paste("Intercept only", sep=""),
                    analysis = list(method="trivial")
                  )
    bf <- list(bf=NA, properror=NA, method=NA)
    BFtry({
      if( nFactors == 0 ){
        numerator = denominator
        bf = list(bf = 0, properror = 0, method = "trivial")
      }else if(all(relevantDataTypes == "continuous")){
        ## Regression
        reg = summary(lm(formula,data=data))
        R2 = reg[[8]]
        N = nrow(data)
        p = length(attr(terms(formula),"term.labels"))
        if( any( names( rscaleEffects ) %in% attr(terms(formula),"term.labels")) ){
          stop("Continuous prior settings set from rscaleEffects; use rscaleCont instead.")
        }
        bf = linearReg.R2stat(N,p,R2,rscale=rscaleCont)
      }else if(all(relevantDataTypes != "continuous")){
        # ANOVA or t test
        freqs <- table(data[[factors[1]]])
        if(all(freqs==1)) stop("not enough observations")
        nLvls <- length(freqs)
        rscale = ifelse(dataTypes[factors[1]] == "fixed", rscaleFixed, rscaleRandom)
        if(length(rscaleEffects)>0)
          if(!is.na(rscaleEffects[factors[1]]))
            rscale = rscaleEffects[factors[1]]
        if( (nFactors==1) & (nLvls==2) ){
          # test
          # independent groups t
          t = t.test(formula = formula,data=data, var.eq=TRUE)$statistic
          bf = ttest.tstat(t=t, n1=freqs[1], n2=freqs[2],rscale=rscale*sqrt(2))
        }else if( (nFactors==1) & (nLvls>2) & all(freqs==freqs[1])){
          # Balanced one-way
          Fstat = summary(aov(formula, data=data))[[1]]["F value"][1,]
          J = length(freqs)
          N = freqs[1]
          bf = oneWayAOV.Fstat(Fstat, N, J, rscale)
        }else if( (nFactors > 1) | ( (nFactors == 1) & any(freqs!=freqs[1]))){
          # Nway ANOVA or unbalanced one-way ANOVA
          bf = nWayFormula(formula=formula, data = data,
                dataTypes = dataTypes,
                rscaleFixed = rscaleFixed,
                rscaleRandom = rscaleRandom,
                rscaleEffects = rscaleEffects,
                posterior = FALSE, ...)
        }else{ # Nothing
          stop("Too few levels in independent variable: ",factors[1])
        }
      }else{
        # GLM
        bf = nWayFormula(formula=formula, data = data,
                       dataTypes = dataTypes,
                       rscaleFixed = rscaleFixed,
                       rscaleRandom = rscaleRandom,
                       rscaleCont = rscaleCont,
                       rscaleEffects = rscaleEffects,
                       posterior = FALSE, ...)
      }
    }) # End try expression

    numerator@analysis = as.list(bf)
    numerator = combineModels(list(numerator,old.numerator))

    bf_df = data.frame(bf = numerator@analysis[['bf']],
                error = numerator@analysis[['properror']],
                time = date(),
                code = randomString(1)
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


