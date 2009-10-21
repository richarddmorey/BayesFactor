printPCLBFAnovaPackage=function(object){
capMod=ifelse(object@settings@Ktype==0,"Cowan","Pashler")
display=paste(
                        "WOMBBAT Analysis\n",
			"-------------------------------------\n",
                        "Analysis Name                  : ",object@settings@analysisname,"\n",
                        "Time and Date                  : ",object@settings@time,"\n",
			"Deviance Information Criterion : ",round(object@output@DIC[1],1)," (Effective parameters: ",round(object@output@DIC[2],1),")\n",
                        "Capacity Model:                : ",capMod,"\n",                       
                        "(view summary for more details)\n",
                        "\n\n",sep="")
cat(noquote(display))
invisible(object)
}

plot.PCLBFAnovaPackage=function(x,...)
{
  plot(x@output@Effchains,...)
}


setClass("PCLBFAnovaSettings", 
      representation(
	data = "data.frame",      # the data to be analyzed, in raw form
	filename = "character",	  # source of the data
	analysisname="character", # name of the analysis, for future reference 
	time="character",         # time and date the analysis was started
	SelCols="list",           # Result of column selection query
	lvlnames="list",          # names of the factor levels of the columns of interest
	newDat="data.frame",      # data, restricted to only columns of interest and collapsed
	mods="list",              # available model effects (before selection)
	effs="data.frame",        # a nice matrix of all possible model effects
	SelEffs="list",		  # result of the model building query
	intMods="data.frame",     # The effects of interest 
	SelEffs2="data.frame",    # nice matrix giving the model requested
	newDat2="data.frame",  # data including the categorical components
	namedDat2="data.frame",   # data with nice names instead of integers for factor levels
	Lvls="integer",             # How many levels do the effects have?
	PriorSetup="list",        # result of prior selection query
	MCMCSetup="list",          # result of MCMC setup query
	designMatrix="matrix"
	),
contains="list")


setClass("PCLBFAnovaOutput", 
      representation(
	chains="mcmc",
	likeChain="numeric",
	par="data.frame",
	burnin="numeric",
	BayesFactor="numeric",
	BayesFactorCI = "numeric",
        predVals="array"                    
),
contains="list")


setClass("PCLBFAnovaPackage",
      representation(
	settings = "list",
	output   = "list"
      ),
contains="list")


setMethod("show", "PCLBFAnovaPackage", printPCLBFAnovaPackage)


summary.PCLBFAnovaPackage=function(object,...){
  capMod=ifelse(object@settings@Ktype==0,"Cowan","Pashler")
  useCovMod = ifelse(object@settings@useCov,"Yes","No")
  time = gsub('\\s','-',object@settings@time,perl=T)
  time = gsub(':','.',time,perl=T)
  name = gsub('[\\s:\\\\]','',object@settings@analysisname,perl=T)

  ## Model stuff
  onK = paste("K        = ",paste(c("muK",as.character(object@settings@SelEffs2[as.logical(object@settings@SelEffs2[,5]),4])),collapse=" + "),collapse="")
  onZ = paste("logit(Z) = ",paste(c("muZ",as.character(object@settings@SelEffs2[as.logical(object@settings@SelEffs2[,6]),4])),collapse=" + "),collapse="")
  onG = paste("logit(G) = ",paste(c("muG",as.character(object@settings@SelEffs2[as.logical(object@settings@SelEffs2[,7]),4])),collapse=" + "),collapse="")

  columns = paste(paste(object@settings@SelCols$selectedcols[,1],object@settings@SelCols$selectedcols[,2],sep=": "),collapse="\n")

  covnames=list()
  if(object@settings@useCov){
    covs = 1:object@settings@nCovMat

    for(i in 1:object@settings@nCovMat){
      x = object@settings@myCovList[lapply(object@settings@CovSetup,function(v) v[[2]])[[i]],c(1,3)]
      covnames[[i]]=paste(x[,1],x[,2],sep=" on ")
      covs[i]=paste(covnames[[i]],collapse=", ")
    }
    covs=paste("  Covariance Matrices:\n",paste(covs,collapse="\n"))
    wishartdf = paste("Wishart df     : ",object@settings@PriorSetup$WishartDF,"\n",sep="")

  }else{
    covs=""
    wishartdf = ""
  }

  if(object@settings@MCMCSetup$useMH=="1"){
  	MCMCtext = paste(
    		   "\nMCMC Settings\n-------------------\n",
    		   "Type           : Random Walk Metropolis-Hastings\n",
		   "Iterations     : ",object@settings@MCMCSetup$nIter,"\n",
    		   "Thin           : ",object@settings@MCMCSetup$MHthin,"\n",
    		   "Effective Iters: ",object@settings@EffectiveIters,"\n",
    		   "Scale          : ",object@settings@MCMCSetup$MHscale,"\n\n",
    sep="")

  }else{
	MCMCtext = paste(
    		   "\nMCMC Settings\n-------------------\n",
    		   "Type           : Hybrid\n",
		   "Iterations     : ",object@settings@MCMCSetup$nIter,"\n",
    		   "Burnin         : ",object@settings@MCMCSetup$burnin,"\n",
    		   "Epsilon        : (",object@settings@MCMCSetup$epsLow,", ",object@settings@MCMCSetup$epsUpp,")\n",
    		   "Leapfrog steps : ",object@settings@MCMCSetup$leapfrog,"\n\n",
    sep="")
  }

  ## Write info
  outputInfo = paste(
    "WOMMBAT Analysis\n",
    "--------------------\n",
    "Analysis name: ",object@settings@analysisname,"\n",
    "Analysis time: ",object@settings@time,"\n",
    "Filename     : ",object@settings@filename,"\n",
    "Model Type   : ",capMod,"\n",
    "Covariances? : ",useCovMod,"\n",
    "\nPrior Specification\n-------------------\n",
    "Inverse Gamma  : (a=",object@settings@PriorSetup$IGa0,", b=",object@settings@PriorSetup$IGb0,")\n",
    "muK            : Normal(mu=",object@settings@PriorSetup$meanMuK,", sigma=",object@settings@PriorSetup$sdMuK,")\n",
    "muZ            : Normal(mu=",object@settings@PriorSetup$meanMuA,", sigma=",object@settings@PriorSetup$sdMuA,")\n",
    "muG            : Normal(mu=",object@settings@PriorSetup$meanMuG,", sigma=",object@settings@PriorSetup$sdMuG,")\n",
    wishartdf, MCMCtext,

    "\nModel\n-------------------\n",
    columns,"\n\n",
    onK,"\n",onZ,"\n",onG,"\n",
    covs,"\n\n",


    "DIC            : ",round(object@output@DIC[1],1)," (Effective parameters: ",round(object@output@DIC[2],1),")\n",
    
    sep="")

cat(noquote(outputInfo))  
invisible(object)
}


