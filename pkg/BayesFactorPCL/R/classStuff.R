


setClass("PCLBF",
	representation(
			t.value="numeric",
			N1="numeric",
			N2="numeric",
			num.samples="numeric",
			JZS="numeric",
			hybrid="numeric",
			NOH="numeric",
			null.d="numeric",
			prior.cauchy="logical",
			prior.scale="numeric",
			pi0="numeric",
			MCMC="logical",
			marginal.d="function",
			MCMC.iters="numeric",
			MCMC.95CI="numeric",
			MCMC.acceptance="numeric",
			MCMC.chains="numeric",
			freq.CI="numeric",
			freq.conf="numeric",
			freq.p="numeric"
	),
contains="list")



PCLBFShow=function(object){
priorType = ifelse(object@prior.cauchy,"Cauchy","Normal")
priorType = paste(priorType," (scale of ",object@prior.scale,") ",sep="")
intType   = ifelse(object@MCMC,paste("MCMC (",object@MCMC.iters," iterations)",sep=""),"Adaptive Quadrature")
mcmc95    = ifelse(object@MCMC,paste("*The MCMC 95% CI on the NOH Bayes Factor is (",zapsmall(object@MCMC.95CI[1]),",",zapsmall(object@MCMC.95CI[2]),")\n",sep=""),"")
display=paste(
			"Bayesian Analysis\n-----------------\n",
			"Prior type         : ",priorType,"\n",
			"JZS Bayes Factor   : ", object@JZS,"\n",
			"NOH Bayes Factor   : ", object@NOH," (null region is ",-object@null.d,",",object@null.d,")\n",
			mcmc95,
			"Hybrid Bayes Factor: ", object@hybrid," (pi0 is ",object@pi0,"; null region is ",-object@null.d,",",object@null.d,")\n",
			"\n\nIntegration Type: ",intType,
			"\n\n",

			"Frequentist Analysis\n--------------------\n",
			100*object@freq.conf,"% CI on d: (",zapsmall(object@freq.CI[1]),",",zapsmall(object@freq.CI[3]),")\n",
			"p value    : ",zapsmall(object@freq.p),			
		
		"\n\n",sep="")
cat(noquote(display))
return(NULL)
}		



setMethod("show", "PCLBF", PCLBFShow)
