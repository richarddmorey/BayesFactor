
doModel = function(X,y,Nks,M=1000,burnin=100,progress=0,MCMCCI=.95,rscale=1)
{

        burnin=1:burnin
	N = dim(X)[1]   # Total Observations
        p = sum(Nks)
        K = length(Nks)

	# PArameter groups
	grps=c(0,inverse.rle(list(values=1:K,lengths=Nks)))


	# reserve space for chains, and init starting values
	start.g = rep(1,K)

	chain.beta = matrix(NA,nrow=dim(X)[2],ncol=M)
	chain.g    = matrix(NA,nrow=K,ncol=M)
	chain.sig2 = 1:M * NA

	chain.beta[,1] = jitter(rep(0,p+1))
	chain.sig2[1] = 1
	chain.g[,1] = start.g

	B = diag(c(0,inverse.rle(list(values=1/start.g,lengths=Nks))))
	XtX = tcrossprod(t(X))

	dens = matrix(NA,K,M-1)



	for(m in 2:M)
	{
		# Draw g
		for(k in 1:K)
		{
			sub.bet = chain.beta[grps==k,m-1]
			chain.g[k,m] = rinvgamma(1, Nks[k]/2 + .5,.5 + .5*sum(sub.bet^2)/chain.sig2[m-1])
		}


        	# draw sigma2
		B = diag(c(0,inverse.rle(list(values=1/chain.g[,m],lengths=Nks))))
        	chain.sig2[m] = rinvgamma(1,N/2+p/2,.5*sum((y-X%*%chain.beta[,m-1])^2)+.5*sum(chain.beta[,m-1]^2*diag(B)))


        	# draw beta
		Sigma = chain.sig2[m] * solve(XtX + B)
		Mu    = Sigma%*%t(X)%*%y/chain.sig2[m]
	
		chain.beta[,m] = mvrnorm(1,Mu,Sigma)

        	# Theorem for conditional distributions
        	# of subvectors of multivariate normal vectors
        	# See Moser (or other linear models text)
        	for(k in 1:K)
		{
			subset = grps==k
          		mu1 = Mu[subset]/sqrt(chain.sig2[m])
          		mu2 = Mu[!subset]/sqrt(chain.sig2[m])
          		known.beta = chain.beta[!subset,m]/sqrt(chain.sig2[m])
          		S11 = Sigma[subset,subset]/chain.sig2[m]
          		S22 = Sigma[!subset,!subset]/chain.sig2[m]
          		S12 = Sigma[subset,!subset]/chain.sig2[m]
          		condMean = mu1 + S12%*%solve(S22)%*%(known.beta-mu2)
			condVariance = S11 - S12%*%solve(S22)%*%t(S12)
          		dens[k,m-1] = dmvnorm(mu1*0,condMean,condVariance)
 		}
	}

	if(K==1){
		MCMCN = effectiveSize(dens[,-burnin])
		vars  = var(apply(dens[,-burnin])

		prior=(1/pi)^Nks


		BF = mean(dens[,-burnin])/prior
		CI = c(BF + qnorm((1-MCMCCI)/2)*sqrt(vars/MCMCN),BF - qnorm((1-MCMCCI)/2)*sqrt(vars/MCMCN))

	}else{
		MCMCN = apply(dens[,-burnin],1,effectiveSize)
		vars  = apply(dens[,-burnin],1,var)

		prior=(1/pi)^Nks


		BF = rowMeans(dens[,-burnin])/prior
		CI = rbind(BF + qnorm((1-MCMCCI)/2)*sqrt(vars/MCMCN),BF - qnorm((1-MCMCCI)/2)*sqrt(vars/MCMCN))
	}
	return(list(chain.beta=chain.beta,chain.sig2=chain.sig2,chain.g=chain.g,density=dens,burnin=burnin,BayesFactor=BF,CI=CI))
}  

