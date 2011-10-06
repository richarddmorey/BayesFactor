
multilevel.Gibbs.AR = function(y,treat,subj,iterations=1000,r.scale=1,betaTheta=5,sdMet=.3, progress=TRUE,return.chains=TRUE)
{

	Nobs = table(subj)
	N = length(Nobs)
	
	
	if(length(treat)!=length(y))
		{
			stop("Invalid condition vector: treat.")
		}
	

	if(length(subj)!=length(y))
		{
			stop("Invalid subject vector: subj.")
		}
	
	
	iterations = as.integer(iterations)

	if(progress){
		progress = round(iterations/100)
		pb = txtProgressBar(min = 0, max = as.integer(iterations), style = 3) 
	}else{ 
		pb=NULL 
	}
	
    pbFun = function(samps){ if(progress) setTxtProgressBar(pb, samps)}

	chains = .Call("RgibbsTwoSampleARmulti", y, N, Nobs, treat, r.scale, betaTheta, iterations, sdMet, progress, pbFun, new.env(), package="BayesFactorPCL")
	
	if(progress) close(pb)
	
	chains = t(matrix(chains,ncol=iterations))
	chains = data.frame(chains)
	colnames(chains) = c(
						paste("mu0",1:N,sep="."),
						paste("beta",1:N,sep="."),
						"muBeta",
						"sig2Beta",
						paste("sig2e",1:N,sep="."),
						"g",	
						"theta"
						)
	
	acc = mean(diff(chains$theta)!=0)
	cat("\n","theta acceptance rate:",acc,"\n")
	
	return(list(chains=mcmc(chains),acc=acc))
}


thetaLogLikeARmulti = function(theta, mu0, beta, sig2, g, y, t, Nobs, betaTheta=5)
{
	N = length(mu0)
	cumobs = cumsum(Nobs)
	maxobs = max(Nobs)

	.Call("RthetaLogLikeARmulti", theta, mu0, beta, sig2, g, y, as.integer(N), t, as.integer(Nobs), as.integer(cumobs), as.integer(maxobs), betaTheta, package="BayesFactorPCL")
	
}